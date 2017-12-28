#!/usr/bin/env python
# create_samplesheet.py - given the name of a run folder, create the sample
# sheet for the run.

import os
import sys
import hcidemux.gnomex
import hcidemux.pipelineparams as params

def EraseCommas(s):
	"""Removes all commas in string s."""
	if s:
		return ''.join(s.split(','))
	else:
		return s

class SampleSheetCreator:
	
	def __init__( self, rundirectory, lanes=[] ):
		self.dirname = rundirectory
		self.id = self.dirname.split('/')[-1]
		self.lanes = lanes[:]

	def log( self, values ):
		sys.stderr.write(" ".join(values) + "\n")

	def WriteSampleSheetRow( self, flowcell_barcode, lane, sampleid, genome, sample_barcode, fname, lname, rnum, ofs ):
		# Regular row. Write it to the file.
		row = [ flowcell_barcode,`lane`,sampleid,EraseCommas(genome) or 'None',sample_barcode or '',EraseCommas(fname+' '+lname),'N','RI','Sandy',rnum ]
		ofs.write(','.join(row)+'\n')

	def CreateSampleSheet(self):
		"""Creates sample sheet for multiplexed flow cell, and writes 
		it to the BaseCalls directory within the run directory. 
		Returns full path of sample sheet file name.
		Changes for version 1.8 of CASAVA:
			Added SampleProject column
			altered other column names
			write file to BaseCalls directory
		4/3/2012: now permitting samples without barcodes. These samples will
		have an Null values for the sample.barcode field, and no barcode
		processing will happen for that lane.
		5/16/2012: if only single sample is in a lane, the barcode for
		that sample will be ignored, and it will be treated as if a single
		non-barcoded sample is there. Doing this because sequence quality
		for single barcodes is poor, and many clients just want all the
		reads anyway.
		11/28/2012: External customer Amplicon Express is submitting libraries
		with 8-base custom barcodes. They submit libraries with any number of
		samples ready for sequencing, but don't want to divulge which samples
		are present. If registered user for Amplicon gets a lane, then the master
		barcode sheet for Amplicon will be subsituted for that lane, so their
		master spreadsheet will get used for every library they submit.
		11/12/2013: Updated to handle paired bar codes.
		03/06/2014: Updated to handle rapid runs. These runs are
		registered in GnomEx as a single lane, but actually have 2
		lanes of data (numbered 1 and 2). GNomEx does not support
		this, so the pipeline is going to fake it by detecting these
		flow cells, and adding a second lane (identical to the first)
		in the sample sheet.
		"""
		# Select lanes from flow cell that are bar coded.
		g = hcidemux.gnomex.GNomExConnection()
		connection = g.GnConnect(params.db_user,params.db_password,asdict=False)
		c = connection.cursor()
		try:
			query = """select flowcell.barcode,
	flowcellchannel.number,
	sample.number,
	genomebuild.genomebuildname,
	sample.barcodesequence,
	appuser.firstname, 
	appuser.lastname,
	request.number,
	sample.barcodesequenceb
from flowcell
	join flowcellchannel on flowcellchannel.idflowcell=flowcell.idflowcell
	join sequencelane on sequencelane.idflowcellchannel = 
		flowcellchannel.idflowcellchannel
	join sample on sequencelane.idsample = sample.idsample
	left outer join genomebuild on sequencelane.idgenomebuildalignto = 
		genomebuild.idgenomebuild
	join request on sequencelane.idrequest = request.idrequest
	join appuser on request.idappuser = appuser.idappuser
where flowcellchannel.filename = '%s'\n""" % self.id
			if self.lanes:
			# Append lane selection to where clause.
				query += "and flowcellchannel.number in (%s)\n" % str(','.join(map(str,self.lanes)))
	
			query += "order by flowcellchannel.number, sample.number;"
			#print query
			c.execute(query)
		except pymssql.OperationError:
			sys.stderr.write( query )
			raise
		results = c.fetchall()
		# Open file
		#ofs = open(samplesheet_fname,'w')
		ofs = sys.stdout
		# Write header.
		header = ['FCID','Lane','SampleID','SampleRef','Index','Description','Control','Recipe','Operator','SampleProject']
		ofs.write(','.join(header)+'\n')
		# Count how many samples are in each lane.
		sample_count = {}
		for rec in results:
			try:
				lane = rec[1]
			except KeyError:
				print rec
				raise
			try:
				sample_count[lane]+=1
			except KeyError:
				sample_count[lane] = 1

		# Write the results.
		for rec in results:
			# Write row.
			flowcell_barcode = rec[0]
			lane = rec[1]
			sampleid=rec[2]
			genome=rec[3]
			sample_barcode=rec[4]
			sample_barcode_b=rec[8]
			if sample_barcode_b is not None:
				sample_barcode = sample_barcode.strip()+'-'+sample_barcode_b.strip()
			elif type(sample_barcode) == type(''):
				sample_barcode = sample_barcode.strip()
			fname=rec[5]
			lname=rec[6]
			rnum=rec[7]
			# If lane has a single sample, don't write its
			# barcode.
			if sample_count[lane] == 1:
				sample_barcode = ''
			self.WriteSampleSheetRow( flowcell_barcode, lane, sampleid, genome, sample_barcode, fname, lname, rnum, ofs )

		# Check if this is a rapid run. This will be the case
		# if the flow cell is registered as a single lane, but there
		# will be two lanes of data for it.
		if len(sample_count.keys())==1 and \
			int(sample_count.keys()[0]) == 1 and \
			os.path.exists(os.path.join(self.dirname,'Data','Intensities','L002')):
			# Its a rapid run.
			self.log(["Run", self.id, "is a rapid run. Duplicating rows for lane 2 in sample sheet",samplesheet_fname,"."])
			for rec in results:
				# Write row.
				flowcell_barcode = rec[0]
				lane = rec[1]
				sampleid=rec[2]
				genome=rec[3]
				sample_barcode=rec[4]
				sample_barcode_b=rec[8]
				if sample_barcode_b is not None:
					sample_barcode = sample_barcode.strip()+'-'+sample_barcode_b.strip()
				elif type(sample_barcode) == type(''):
					sample_barcode = sample_barcode.strip()
				fname=rec[5]
				lname=rec[6]
				rnum=rec[7]
				# If lane has a single sample, don't write its
				# barcode.
				if sample_count[lane] == 1:
					sample_barcode = ''
				self.WriteSampleSheetRow( flowcell_barcode, 2, sampleid, genome, sample_barcode, fname, lname, rnum, ofs )



		# Close file.
		#ofs.close()

		# Close database.
		connection.close()
		# Return file name.

if __name__ == "__main__":
	if len(sys.argv) < 2 or sys.argv[1] == '--help':
		sys.stderr.write("%s creates a sample sheet for an Illumina\n" % sys.argv[0] )
		sys.stderr.write("run folder from information in GNomEx and writes it to stdout.\n" )
		sys.stderr.write("Use: %s <run_folder> [lane# lane# ... ]\n" % sys.argv[0] )
		sys.exit(1)
	s = SampleSheetCreator( sys.argv[1],map(int,sys.argv[2:]))
	s.CreateSampleSheet()
