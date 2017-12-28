import os
import re
import sys
import subprocess

import fastq
import pipelineparams as params
from states import States
from run import Run

class PatchPcrRun(Run):

	def __init__( self, id, full_path_of_run_dir, state=None ):
		Run.__init__(self, id, full_path_of_run_dir, state )
		self.type_of_runs = "PatchPcr"
		self.transition = {
			# curr_state, function, success_state, fail_state
			States.new: (PatchPcrRun.CheckRegisteredVerbose,
				States.wait_for_data,
				States.check_registered ),
			States.check_registered: (PatchPcrRun.CheckRegisteredSilent,
				States.wait_for_data,
				States.check_registered),
			States.wait_for_data: (PatchPcrRun.CheckTransferComplete,
				States.check_if_miseq,
				States.wait_for_data),
			States.check_if_miseq: (PatchPcrRun.CheckIfMiseq,
				States.check_if_patchpcr,
				States.error_detected),
			States.check_if_patchpcr: (PatchPcrRun.CheckIfPatchPcr,
				States.patchpcr_sample_sheet,
				States.error_detected),
			States.patchpcr_sample_sheet: (PatchPcrRun.PatchPcrSampleSheet,
				States.patchpcr_demultiplex,
				States.error_detected),
			States.patchpcr_demultiplex: (PatchPcrRun.PatchPcrDemultiplex,
				States.patchpcr_postprocess,
				States.error_detected),
			States.patchpcr_postprocess: (PatchPcrRun.PatchPcrPostprocess,
				States.patchpcr_distribute,
				States.error_detected),
			States.patchpcr_distribute: (PatchPcrRun.PatchPcrDistribute,
				States.Miseq_qc,
				States.error_detected),
			States.qc: (PatchPcrRun.Qc_18,
				States.archive,
				States.error_detected),
			States.archive: (PatchPcrRun.Archive,
				States.complete,
				States.error_detected),
			States.error_detected: (PatchPcrRun.NotifyError,
				States.error_notified,
				States.error_detected),
		}

	def PatchPcrDistribute( self ):
		"""
		Copies data from Patch PCR runs to GNomEx.
		"""
		self.Log("Distributing files for run %s" % self.id)
		self.CopyDataFiles()
		return True

	def CheckIfPatchPcr( self ):
		"""
		Determines if a miseq run is a patch pcr run.
		"""
		#self.Log("Checking if run %s is a patch pcr run." % self.id)
		query = """select distinct application.codeapplication, flowcellchannel.filename, application.application
from application 
join seqlibprotocolapplication 
	on application.codeapplication = seqlibprotocolapplication.codeapplication
join seqlibprotocol on seqlibprotocol.idseqlibprotocol = seqlibprotocolapplication.idseqlibprotocol
join sample on sample.idseqlibprotocol = seqlibprotocol.idseqlibprotocol
join sequencelane on sequencelane.idsample = sample.idsample
join flowcellchannel on sequencelane.idflowcellchannel = flowcellchannel.idflowcellchannel
where application.application like '%%Patch PCR%%' and 
	flowcellchannel.filename = '%s';""" % self.id
		c = self.Query(query)
		results = c.fetchall()
		is_patchpcr = len(results)>0
		if is_patchpcr:
			self.Log("Run %s is a patch pcr run." % self.id)
		#else:
		#	self.Log("Run %s is not a patch pcr run." % self.id)
		return is_patchpcr



	def PatchPcrSampleSheet( self ):
		"""Creates sample sheet for multiplexed flow cell, and writes 
		it to the BaseCalls directory within the run directory. 
		Returns full path of sample sheet file name.
		"""
		self.Log("Generating patch PCR run sample sheet for run %s." % self.id )
		# Generate sample sheet name. If it already exists just
		# return its name. Doing this to make it easier to override
		# the automatic sample sheet creation. Brett Milash, 10/11/2013.
		self.sample_sheet = self.SampleSheetName()
		if os.path.exists(self.sample_sheet):
			self.Log("Sample sheet for run %s already exists." % self.id)
			return True

		# Select lanes from flow cell that are bar coded.
		try:
			query = """select sample.number,
	sample.barcodesequence,
	request.number
from flowcell
	join flowcellchannel on flowcellchannel.idflowcell=flowcell.idflowcell
	join sequencelane on sequencelane.idflowcellchannel = 
		flowcellchannel.idflowcellchannel
	join sample on sequencelane.idsample = sample.idsample
	left outer join genomebuild on sequencelane.idgenomebuildalignto = 
		genomebuild.idgenomebuild
	join request on sequencelane.idrequest = request.idrequest
	join appuser on request.idappuser = appuser.idappuser
where flowcellchannel.filename = '%s'
order by flowcellchannel.number, sample.number;""" % self.id
			c = self.Query( query )
		except pymssql.OperationError:
			print query
			raise
		results = c.fetchall()
		# Open file
		ofs = open(self.sample_sheet,'w')
		# Copy the header from the existing sample sheet.
		ifs = open(os.path.join(self.dirname,'SampleSheet.csv'))
		rec = ifs.readline()
		while rec and not rec.startswith('Sample_ID'):
			ofs.write(rec)
			rec = ifs.readline()
		# Write the Sample_ID header.
		ofs.write(rec)
		ifs.close()

		# Write the results.
		s = [''] * 10
		for rec in results:
			# Write row.
			# Sample id.
			s[0]=rec[0]
			# Barcode.
			s[5]=rec[1]
			# Request number
			s[8]=rec[2]
			ofs.write(','.join(s)+ '\n')

		# Close file.
		ofs.close()
		return True


	def PatchPcrDemultiplex( self ):
		"""
		Calls bcl2fastq2 program to demultiplex the run.
		"""
		self.Log("Demultiplexing patch PCR run %s." % self.id)
		# Set minimum read length to the length of the bar code 
		# read carrying the random n-mer. This will always be the 
		# 3rd read.
		read_lengths = self.ReadLengths()
		min_read_length = read_lengths[2]
		if self.CheckIfPatchPcr():
		# Determine the use-bases-mask.
			if len(read_lengths) == 4:
			# This is a paired-end run with two index reads.
				use_bases_mask = "Y*,I*,Y*,Y*"
			else:
			# This is a single-end run with two index reads.
				use_bases_mask = "Y*,I*,Y*"



		child_dir=self.dirname
		output_dir = os.path.join( self.dirname, "Unaligned" )
		child_stdout=open(os.path.join(child_dir,"bcl2fastq.out"),'w')
		child_stderr=open(os.path.join(child_dir,"bcl2fastq.err"),'w')
		child_args=[os.path.join(params.bcl2fastq2_dir,"bcl2fastq"),
			"--runfolder-dir",self.dirname,
			"--output-dir",output_dir,
			"--sample-sheet",self.sample_sheet,
			"--barcode-mismatches","1",
			"--use-bases-mask",use_bases_mask,
			"--minimum-trimmed-read-length", `min_read_length`
		]
		p = subprocess.Popen(args=child_args,cwd=child_dir,stdout=child_stdout,stderr=child_stderr)
		retval = p.wait()
		child_stdout.close()
		child_stderr.close()
		return True

	def PatchPcrPostProcessSample(self,project,sample):
		"""
		Postprocesses data files for one sample in a Patch PCR run.
		Takes the sequence from the read 2 file and appends it to
		the read names in the read 1 and read 3 files (if present).
		"""
		self.Log("Post processing sample %s." % sample)
		# Create input file names for reads R1, R2, and R3.
		r1_in_file = os.path.join(self.dirname,"Unaligned",project,"%s_L001_R1_001.fastq.gz"%sample)
		r2_in_file = os.path.join(self.dirname,"Unaligned",project,"%s_L001_R2_001.fastq.gz"%sample)
		r3_in_file = os.path.join(self.dirname,"Unaligned",project,"%s_L001_R3_001.fastq.gz"%sample)

		# Create output file names for R1 and R2.
		gnomex_sample=sample.split("_")[0]
		r1_out_file = os.path.join(self.dirname,"Unaligned","%s_%s_1_1.txt.gz" % ( gnomex_sample, self.id))
		r2_out_file = os.path.join(self.dirname,"Unaligned","%s_%s_1_2.txt.gz" % ( gnomex_sample, self.id))

		# Process read 1.
		r1_ifs=fastq.Reader(r1_in_file)
		r2_ifs=fastq.Reader(r2_in_file)
		r1_ofs=fastq.Writer(r1_out_file)
		for r1_read in r1_ifs:
			nmer_read = r2_ifs.next()
			r1_read.name += "-" + nmer_read.sequence
			r1_ofs.write(r1_read)
		r1_ifs.close()
		r2_ifs.close()
		r1_ofs.close()
		self.datafiles.append(r1_out_file)

		# If read 3 exists, process it.
		if os.path.exists(r3_in_file):
			# Process read 1.
			r3_ifs=fastq.Reader(r3_in_file)
			r2_ifs=fastq.Reader(r2_in_file)
			r2_ofs=fastq.Writer(r2_out_file)
			for r3_read in r3_ifs:
				nmer_read = r2_ifs.next()
				r3_read.name += "-" + nmer_read.sequence
				r2_ofs.write(r3_read)
			r3_ifs.close()
			r2_ifs.close()
			r2_ofs.close()
			self.datafiles.append(r2_out_file)

	def PatchPcrPostprocess( self ):
		"""
		Postprocessing for patch PCR runs. Takes the random n-mer, which
		is produced by the second barcode read file, and should be in
		the file with the _2 suffix, and appends the random n-mer sequence
		to the read names of the data read(s).
		"""
		# Get projects and samples for this run from the sample sheet.
		projects=self.GetSampleSheetProjectsSamples()

		# For each project...
		for project in projects.keys():
			# For each sample...
			for sample in projects[project]:
				self.PatchPcrPostProcessSample(project,sample)

		# Generate the MD5 checksums.
		return self.GenerateChecksums()

def main():
	p_run = PatchPcrRun("bar","/foo/bar")

if __name__ == "__main__":
	main()
