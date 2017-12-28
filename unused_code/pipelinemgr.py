#!/usr/bin/python
"""PipelineMgr.py - program to manage running the HiSeq data processing pipeline.
Consists of a PipelineMgr class with methods to move runs through the pipeline.
Each step in pipeline is defined by a state in a state machine. A run's current
state defines which method to use, and what the possible next states could be.

2013/06/11 - replaced /bin/cp with /usr/bin/rsync to prevent data corruption
	when copying files to data repository.

"""

try:
	import sqlite3		# For object persistence.
except ImportError:
	import sqlite as sqlite3
import xml.dom.minidom
import time
import os
import socket
import sys
import re
import string
import GNomEx
import subprocess
import traceback
from Emailer import Emailer
from Logger import Logger
from SimultaneousJobRunner import SimultaneousJobRunner
import Fastq

sys.path.append(os.getcwd())
import AmpliconExpress
import PipelineParams as params

# Connect to the datababase. This should be global, so it won't get 
# garbarge collected when a function ends.
connection = GNomEx.GNomExConnection().GnConnect(params.db_user,params.db_password)

def EraseCommas(s):
	"""Removes all commas in string s."""
	if s:
		return ''.join(s.split(','))
	else:
		return s

# Import the module for System V semaphores.
try:
	import sysv_ipc
except ImportError:
	Logger().Log( "%s requires module sysv_ipc. See http://semanchuk.com/philip/sysv_ipc/." % sys.argv[0] )
	sys.exit(1)

class States:
	"""States: legal states in the state machine. The state values are
	strings so they can be displayed meaningfully (as opposed to numbers)
	and can be ordered by their values."""
	new='a_new'
	check_registered='b_check_registered'
	wait_for_data='c_wait_for_data'
	make_sample_sheet='d_make_sample_sheet'
	check_single_index_length='e_check_single_index_length'
	simple_bcl_convert='f_simple_bcl_convert'
	complex_bcl_convert='g_complex_bcl_convert'
	merge_bcl_conversions='h_merge_bcl_conversions'
	distribute_demultiplexed='j_distribute_demultiplexed'
	cleanup ='k_cleanup' 
	qc='l_qc'
	archive='m_archive'
	complete='n_complete'
	error_detected='o_error_detected'
	error_notified='p_error_notified'
	reprocess='q_reprocess'
	check_if_miseq='r_check_if_miseq'
	check_if_patchpcr='s_check_if_patchpcr'
	patchpcr_sample_sheet='t_patchpcr_sample_sheet'
	patchpcr_demultiplex='u_patchpcr_demultiplex'
	patchpcr_postprocess='v_patchpcr_postprocess'
	patchpcr_distribute='w_patchpcr_distribute'

class Run(Logger,Emailer,SimultaneousJobRunner):
	"""Information about one sequencing run including its current
	state in the pipeline."""

	def __init__( self, id, full_path_of_run_dir, state=None ):
		self.id = id
		self.dirname = full_path_of_run_dir
		if state is None:
			self.state = States.new
		else:
			self.state = state
		self.runtype = None
		self.sample_sheet = None
		# List of full pathnames of all data files produced by this run, following rename.
		self.datafiles=[]
		self.output_dirs=[]
		# List of Lane objects.
		self.lanes = []
		SimultaneousJobRunner.__init__( self, verbose=True )
		self.InitializeCoreFacility()
	
	def __cmp__( self, other ):
		"""Sort runs by decreasing state, increasing run start date."""
		if self.state > other.state:
			return -1
		elif self.state < other.state:
			return 1
		elif self.id < other.id:
			return -1
		else:
			return 1

	def Query( self, query ):
		"""
		Runs a SQL query against GNomEx, and retuns a cursor from
		which the results can be retrieved.
		"""
		global connection
		c = connection.cursor()
		c.execute(query)
		return c
	
	def InitializeCoreFacility( self ):
		"""
		InitializeCoreFacility identifies which core facility
		started this run.
		"""
		self.corefacilityname = None
		query="""select distinct facilityname
		from corefacility
		join flowcell on flowcell.idcorefacility = corefacility.idcorefacility
		join flowcellchannel on flowcellchannel.idflowcell = flowcell.idflowcell
		and flowcellchannel.filename = '%s';""" % self.id
		recs=self.Query(query).fetchall()
		try:
			self.corefacilityname = recs[0][0]
		except:
			self.corefacilityname = 'Unknown'
		self.Log("Run %s came from core facility '%s'." % (self.id,self.corefacilityname))

	def FromAddr( self ):
		"""
		Returns from address for emailing purposes.
		"""
		# Unix account under which PipelineMgr.py runs.
		username = os.environ['LOGNAME']
		return username + "@" + socket.getfqdn()

	def Notify( self, subject, message ):
		to_addr = params.notify_addresses
		self.send_msg( self.FromAddr(), to_addr, subject, message, attachments=[] )

	def NotifyError( self ):
		"""Notify someone that this run has entered the error state."""
		to_addr = params.notify_addresses
		subject = "HiSeq Pipeline Notification"
		message = "Run %s (%s) just entered error state." % (self.id,self.dirname)
		self.Log(["Notifying",to_addr,":",message])
		try:
			self.send_msg( self.FromAddr(), to_addr, subject, message, attachments=[] )
			return True
		except:
			self.Log(["PROBLEM! Can't notify",to_addr,"that",message])
			return False

	def CheckRegisteredSilent( self ):
		"""Boolean function to test if run folder actually registered in GNomEx."""
		g = GNomEx.GNomExConnection()
		connection = g.GnConnect(params.db_user,params.db_password)
		c = connection.cursor()
		c.execute("select count(*) from flowcellchannel where filename = %s", self.id)
		results = c.fetchall()
		numlanes = results[0][0]
		# 2013/07/15: Flow cells can now either have 8 lanes (HiSeq2000)
		# or 1 lane (HiSeq2500). Permit either of these styles of
		# flow cell here. Brett Milash.
		if numlanes == 8 or numlanes == 1:
			return True
		else:
			# 2013/10/02 - If run is complete and flow cell still
			# not registered, send an urgent email.
			if self.CheckTransferComplete(notify=False):
				subject = "HiSeq Pipeline Problem! (%s)" % self.id
				message = "The sequencing run %s is done, but the run folder has not been entered in GNomEx.\n" % self.id
				message += "The pipeline can not start until the flow cell has been entered into GNomEx.\n"
				self.NotifyLabStaff( subject, message )
			return False

	def NotifyLabStaff( self, subject, message ):
		# Notify the lab folks to correct the problem.
		to_addr = params.lab_staff_addresses[self.corefacilityname]
		self.send_msg( self.FromAddr(), to_addr, subject, message )

	def CheckRegisteredVerbose( self ):
		if self.CheckRegisteredSilent():
			return True
		else:
			# Complain to lab staff, and log it.
			self.Log(["PROBLEM! Run directory",self.id,"not registered correctly in GNomEx."])
			# Notify the lab folks to correct the problem.
			subject = "HiSeq Pipeline Problem! (%s)" % self.id
			message = "The HiSeq run folder %s has not been entered in GNomEx.\n" % self.id
			message += "Please check GNomEx to make sure the run folder was entered correctly.\n"
			message += "Thanks, and have a nice day!\n"
			self.NotifyLabStaff( subject, message )
			return False

	def CheckTransferComplete(self,notify=True):
		"""Checks if data transfer from sequencer is done. Illumina 
		RTA 1.12 creates a data transfer complete file for each read
		in the run. This function checks if the last file for the
		run is present."""
		configfile = os.path.join(self.dirname, "RunInfo.xml")
		doc = xml.dom.minidom.parse(configfile)
		numreads=len(doc.getElementsByTagName("Read"))
		done_file=os.path.join(self.dirname,"Basecalling_Netcopy_complete_Read%d.txt"%numreads)
		if os.path.exists(done_file):
			self.Log(["Run", self.id, "complete-file",done_file,"exists."])
			if notify:
				self.NotifyStarting()
			return True
		else:
			self.Log(["Run", self.id, "complete-file",done_file,"doesn't exist."])
			return False
	
	def CheckSingleIndexLength( self ):
		if self.sample_sheet is None:
			# Try to find the sample sheet. It should exist.
			self.sample_sheet = self.SampleSheetName()

		# Open, read, and parse the sample sheet. The Index column lists the bar codes. Determine if they are all
		# the same length or not. If they are, return True, else False.
		ifs = open(self.sample_sheet)
		lines = ifs.readlines()
		ifs.close()
		lengths=set()
		for line in lines[1:]:
			f = line.split(',')
			barcode = f[4]
			lengths.add(len(barcode))
		# Count the number of elements in set lengths.
		numlengths=0
		while lengths:
			lengths.pop()
			numlengths+=1
		if numlengths == 1:
			return True
		return False

	def MakeSampleSheet( self ):
		self.sample_sheet = self.CreateSampleSheet()
		if self.sample_sheet is None:
			self.Log(["PROBLEM! Can't locate or create sample sheet .csv file in", self.dirname, "for barcoded run", self.id] )
			return False
		return True
	
	def ReadLengths( self ):
		"""
		Returns a tuple with the lengths of each read in the run 
		read from the RunInfo.xml file.
		"""
		lengths=[]
		doc = xml.dom.minidom.parse(os.path.join(self.dirname,"RunInfo.xml"))
		reads=doc.getElementsByTagName("Read")
		for read in reads:
			lengths.append(int(read.getAttribute("NumCycles")))
		return tuple(lengths)

	def FlowCellVersion( self ):
		"""Returns run's flow cell version from the 
		runParameters.xml file."""
		# Look in the runParameters.xml document for the Flowcell element to determine flow cell type.
		doc = xml.dom.minidom.parse(os.path.join(self.dirname,"runParameters.xml"))
		e=doc.getElementsByTagName("Flowcell")[0]
		flowcell=e.firstChild.nodeValue
		return flowcell

	def BclConvert( self, sample_sheet=None, output_dir=None, use_bases_mask=None, compress_bcls=True, mismatches=1 ):
		"""Converts the bcl files into compressed Fastq."""
		# Generate a sample sheet if it is not already present.
		if sample_sheet:
			self.sample_sheet = sample_sheet
		elif not self.sample_sheet:
			self.sample_sheet = self.SampleSheetName()
		if self.sample_sheet is None:
			self.Log(["PROBLEM! Can't locate or create sample sheet .csv file in", self.dirname, "for barcoded run", self.id] )
			return False

		flowcell = self.FlowCellVersion()
		if flowcell == "HiSeq Flow Cell v3":
			# V3 flowcell processing.
			# Run the configureBclToFastq.pl scripts.
			self.Log(["Running configureBclToFastq on run", self.id])
			child_dir=self.dirname
			if not output_dir:
				output_dir = os.path.join( self.dirname, "Unaligned" )
			child_stdout=open(os.path.join(child_dir,"bcltofastq.out"),'w')
			child_stderr=open(os.path.join(child_dir,"bcltofastq.err"),'w')
			child_args=[os.path.join(params.casava_dir,"configureBclToFastq.pl"),
				"--input-dir",os.path.join(self.dirname,"Data","Intensities","BaseCalls"),
				"--positions-format",".clocs",
				"--sample-sheet",self.sample_sheet,
				"--output-dir",output_dir,
				"--mismatches",`mismatches`,
			]
			if use_bases_mask:
				child_args.append( "--use-bases-mask" )
				child_args.append( use_bases_mask )
			p = subprocess.Popen(args=child_args,cwd=child_dir,stdout=child_stdout,stderr=child_stderr)
			retval = p.wait()
			child_stdout.close()
			child_stderr.close()
		else:
			# V4 and other types of flowcell processing.
			# Run the configureBclToFastq.pl scripts.
			self.Log(["Running configureBclToFastq on run", self.id])
			#basecallsdir=$rundir/Data/Intensities/BaseCalls
			#configureBclToFastq.pl --input-dir=$basecallsdir \
			#	--output-dir=$rundir/Unaligned \
			#	--sample-sheet=$basecallsdir/created_samplesheet.csv \
			#	--positions-format=.clocs \
			#	--mismatches 1 \
			#	--no-eamss
			child_dir=self.dirname
			if not output_dir:
				output_dir = os.path.join( self.dirname, "Unaligned" )
			child_stdout=open(os.path.join(child_dir,"bcltofastq.out"),'w')
			child_stderr=open(os.path.join(child_dir,"bcltofastq.err"),'w')
			child_args=[os.path.join(params.bcl2fastq_dir,"configureBclToFastq.pl"),
				"--input-dir",os.path.join(self.dirname,"Data","Intensities","BaseCalls"),
				"--positions-format",".clocs",
				"--sample-sheet",self.sample_sheet,
				"--output-dir",output_dir,
				"--mismatches",`mismatches`,
				"--ignore-missing-bcl",
				"--ignore-missing-stats",
				"--no-eamss",
			]
			if use_bases_mask:
				child_args.append( "--use-bases-mask" )
				child_args.append( use_bases_mask )
			p = subprocess.Popen(args=child_args,cwd=child_dir,stdout=child_stdout,stderr=child_stderr)
			retval = p.wait()
			child_stdout.close()
			child_stderr.close()

		# If configure successful...
		if retval == 0:
			# Run make in the Unaligned directory.
			self.Log(["Running make (to convert bcl files) on run",self.id])
			child_dir = output_dir
			child_stdout=open(os.path.join(child_dir,"make.out"),'w')
			child_stderr=open(os.path.join(child_dir,"make.err"),'w')
			child_args=["/usr/bin/make","-j","8"]
			p = subprocess.Popen(args=child_args,cwd=child_dir,stdout=child_stdout,stderr=child_stderr)
			retval = p.wait()
			child_stdout.close()
			child_stderr.close()
		if retval == 0 and compress_bcls == True:
			# Compress .bcl files.
			self.CompressFiles("*.bcl")
		return retval == 0
	
	def CheckMultiplex( self ):
		"""Checks if a run contains any barcoded samples or not. Counts
		the number of samples for each lane on the flow cell, and returns
		true if any lane has more than one sample."""
		g = GNomEx.GNomExConnection()
		connection = g.GnConnect(params.db_user,params.db_password,asdict=True)
		c = connection.cursor()
		try:
			query="""select count(*) n
				from flowcellchannel
				join sequencelane on sequencelane.idflowcellchannel = flowcellchannel.idflowcellchannel
				join sample on sample.idsample = sequencelane.idsample
				where flowcellchannel.filename = '%s' 
				group by flowcellchannel.filename, flowcellchannel.number""" % self.id
			c.execute(query)
		except pymssql.OperationError:
			print query
			raise

		results = c.fetchall()
		connection.close()
		for rec in results:
			if rec['n'] > 1:
				self.Log("Run %s is barcoded."%self.id)
				return True
		self.Log("Run %s is not barcoded."%self.id)
		return False

	def FindBarcodedLanes( self ):
		"""Determines which lanes on the flow cell contain barcoded samples. Returns a string
		with the lane numbers, for example "12348"."""
		g = GNomEx.GNomExConnection()
		connection = g.GnConnect(params.db_user,params.db_password)
		c = connection.cursor()
		try:
			query="""select flowcellchannel.number, count(*) n
			from flowcellchannel
			join sequencelane on sequencelane.idflowcellchannel = flowcellchannel.idflowcellchannel
			join sample on sample.idsample = sequencelane.idsample
			where flowcellchannel.filename = '%s'
			group by flowcellchannel.filename, flowcellchannel.number   
			order by flowcellchannel.filename, flowcellchannel.number"""%self.id

			c.execute(query)
		except pymssql.OperationError:
			print query
			raise

		results = c.fetchall()
		connection.close()
		barcoded_lanes = ''
		for rec in results:
			lane = rec[0]
			numsamples=rec[1]
			if numsamples > 1:
				barcoded_lanes += str(lane)
		return barcoded_lanes

	def LocateSampleSheet( self ):
		"""Finds sample sheet file, a comma-delimited file (.csv) in the run directory. Returns
		None if not found."""
		fnames = os.listdir(self.dirname)
		for file in os.listdir(self.dirname):
			if file[-4:] == ".csv":
				return os.path.join(self.dirname,file)
		return None

	def WriteSampleSheetRow( self, flowcell_barcode, lane, sampleid, genome, sample_barcode, fname, lname, rnum, ofs ):
		if AmpliconExpress.custom_handling and fname == AmpliconExpress.firstname and lname == AmpliconExpress.lastname:
			# This is a lane from Amplicon Express. Substitute their
			# master samplesheet, and write it to the file.
			sample_sheet = AmpliconExpress.SampleSheet
			# Substitute the flow cell, lane, and request number into the sample sheet.
			sample_sheet=re.sub('<fcid>',flowcell_barcode,sample_sheet)
			sample_sheet=re.sub('<lane>',`lane`,sample_sheet)
			sample_sheet=re.sub('<request>',rnum,sample_sheet)
			ofs.write( sample_sheet )
		else:
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
		# Generate sample sheet name. If it already exists just
		# return its name. Doing this to make it easier to override
		# the automatic sample sheet creation. Brett Milash, 10/11/2013.
		samplesheet_fname = self.SampleSheetName()
		if os.path.exists(samplesheet_fname):
			return samplesheet_fname

		# Select lanes from flow cell that are bar coded.
		g = GNomEx.GNomExConnection()
		connection = g.GnConnect(params.db_user,params.db_password,asdict=True)
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
where flowcellchannel.filename = '%s'
order by flowcellchannel.number, sample.number;""" % self.id
			c.execute(query)
		except pymssql.OperationError:
			print query
			raise
		results = c.fetchall()
		# Open file
		ofs = open(samplesheet_fname,'w')
		# Write header.
		header = ['FCID','Lane','SampleID','SampleRef','Index','Description','Control','Recipe','Operator','SampleProject']
		ofs.write(','.join(header)+'\n')
		# Count how many samples are in each lane.
		sample_count = {}
		for rec in results:
			lane = rec[1]
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
			self.Log(["Run", self.id, "is a rapid run. Duplicating rows for lane 2 in sample sheet",samplesheet_fname,"."])
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
		ofs.close()
		# Close database.
		connection.close()
		# Return file name.
		return samplesheet_fname

	def DeMultiplex( self ):
		"""Performs barcode processing."""
		# Locate sample sheet file, a comma-delimited file (.csv) in the run directory.
		sample_sheet = self.LocateSampleSheet()
		if sample_sheet is None:
			sample_sheet = self.CreateSampleSheet()
		if sample_sheet is None:
			self.Log(["PROBLEM! Can't locate or create sample sheet .csv file in", self.dirname, "for barcoded run", self.id] )
			return False
		# Identify which lanes contain barcoded samples.
		barcoded_lanes = self.FindBarcodedLanes()
		# Create config template file.
		config_template = os.path.join(self.dirname,"config.Template.txt")
		if not os.path.exists( config_template ):
			if not self.CreateGeraldConfig(config_template):
				return False
 		# Run demultiplex.pl script.
		self.Log(["Running demultiplex.pl on run", self.id])
		basecall_dir=os.path.join(self.dirname,"Data","Intensities","BaseCalls")
		demultiplexed_dir=os.path.join(basecall_dir,"Demultiplexed")
		child_dir=os.path.join(self.dirname,"Data")
		child_stdout=open(os.path.join(child_dir,"demultiplex.out"),'w')
		child_stderr=open(os.path.join(child_dir,"demultiplex.err"),'w')
		child_args=[ os.path.join(params.casava_dir,"demultiplex.pl"),
			"--input-dir",basecall_dir,
			"--sample-sheet",sample_sheet,
			"--alignment-config",config_template,
			"--output-dir",demultiplexed_dir,
			"--mismatches", "1",
		]
		p = subprocess.Popen(args=child_args,cwd=child_dir,stdout=child_stdout,stderr=child_stderr)
		retval = p.wait()
		child_stdout.close()
		child_stderr.close()
		
		if retval == 0:
			# Run make in the Demultiplexed directory.
			self.Log(["Running make (to demultiplex qseq files) on run",self.id])
			child_dir = demultiplexed_dir
			child_stdout=open(os.path.join(child_dir,"make.out"),'w')
			child_stderr=open(os.path.join(child_dir,"make.err"),'w')
			child_args=["/usr/bin/make","-j","8","ALIGN=yes"]
			p = subprocess.Popen(args=child_args,cwd=child_dir,stdout=child_stdout,stderr=child_stderr)
			retval = p.wait()
			child_stdout.close()
			child_stderr.close()
		return retval == 0
	
	def FindGeraldDirectory(self,parent_dir=None):
		"""Returns name of GERALD directory for run, or None if
		no directory exists."""
		if parent_dir is None:
			parent_dir = os.path.join(self.dirname,"Data","Intensities","BaseCalls")
		subdirs = os.listdir(parent_dir)
		subdirs.sort()
		subdirs.reverse()
		pattern = re.compile("GERALD_[0-9]*-[0-9]*-[0-9]*_"+params.username)
		for subdir in subdirs:
			if pattern.match(subdir):
				return os.path.join(parent_dir,subdir)
		return None

	def DetermineRunType(self):
		"""
		DetermineRunType determines the run type of the run by looking
		at the RunInfo.xml file in the run directory. This mechanism
		should be more reliable than checking in GNomEx, since the 
		GNomEx database run type names change from time to time.
		Looking at the RunInfo.xml file, this function counts the
		number of reads >= 25 bases in length. If the number of reads
		>= 25 bp is 2, the run is paired, otherwise single-end.
		Brett Milash, 10/09/2013.
		Updated 03/17/2014: Also determines how many index reads
		are present, and stores the number of data reads and
		index reads in member variables num_data_reads and
		num_index_reads. Brett Milash.
		"""
		# Parse the RunInfo.xml file.
		configfile = os.path.join(self.dirname, "RunInfo.xml")
		doc = xml.dom.minidom.parse(configfile)

		# Get the Read elements from the document.
		reads=doc.getElementsByTagName("Read")

		# Count the reads that are at least 25 bp long.
		self.num_data_reads = 0
		self.num_index_reads = 0
		for read in reads:
			if int(read.getAttribute('NumCycles')) >= 25:
				self.num_data_reads += 1
			else:
				self.num_index_reads += 1

		if self.num_data_reads == 1:
			self.IsPaired = False
		elif self.num_data_reads == 2:
			self.IsPaired = True
		else:
			raise Exception("Expected 1 or 2 data reads in RunInfo.xml. Found %d data reads." % self.num_data_reads )


	def CreateGeraldConfig(self,gerald_config_file,lanes="12345678"):
		"""Creates GERALD configuration file."""
		self.DetermineRunType()
		# Determine if single-end or paired-end sequencing.
		if self.runtype == 'Paired-end reads':
			analysis = 'sequence_pair'
		elif self.runtype == 'Single-end reads':
			analysis = 'sequence'
		else:
			# Unknown type of sequencing.
			self.Log(["PROBLEM! Run ",self.id,"has unknown type of sequencing,",self.runtype,", unable to create GERALD configuration file."])
			return False

		ofs = open(gerald_config_file,'w')
		ofs.write("""# %s
# Created by %s.
%s:ANALYSIS %s
""" % ( gerald_config_file, sys.argv[0], lanes, analysis ) )
		# This section not necessary. Can just use a generic
		# gerald config file to process all lanes, barcoded or not.
		## Write additional lines if there are any non-barcoded lanes.
		#if len(lanes) != 8:
		#	# Identify non-barcoded lanes.
		#	non_barcoded_lanes = ''
		#	for lane in [ '1','2','3','4','5','6','7','8',]:
		#		if lanes.find(lane) == -1:
		#			non_barcoded_lanes += lane
		#	ofs.write("SAMPLE unknown %s:ANALYSIS %s\n" % ( non_barcoded_lanes, analysis ))
		ofs.write("""USE_BASES all
EMAIL_LIST brett.milash@hci.utah.edu
EMAIL_DOMAIN hci.utah.edu
EMAIL_SERVER hci-mail.hci.utah.edu:25
WITH_SEQUENCE true
""")
		ofs.close()
		return True

	def FindUnbarcodedLanes( self, sample_sheet_file ):
		"""Reads sample sheet file to identify lanes that do NOT include
		barcoded samples. Return these lanes as a list of integers."""
		barcoded_lanes = []
		ifs = open(sample_sheet_file)
		for rec in ifs:
			try:
				lane = int(rec.strip().split(',')[1])
				barcoded_lanes.append( lane )
			except ValueError:
				pass
		ifs.close()
		# barcoded_lanes is list of lanes with barcoded samples.
		s1 = set(range(1,9))
		b = set(barcoded_lanes)
		non_barcoded_set = s1 - b
		return list(non_barcoded_set)

	def FindSampleSheetFile( self ):
		for filename in os.listdir(self.dirname):
			if filename[-4:] == ".csv":
				return os.path.join(self.dirname,filename)
		raise IOError, "Sample sheet not found in " + self.dirname

	def MakeBinList( self ):
		"""Finds list of subdirectories named 001, 002, ...."""
		binlist = []
		for entry in os.listdir("."):
			if os.path.isdir(entry) and entry[0] == "0":
				binlist.append(entry)
		binlist.sort()
		return binlist

	def ProcessSample(self,lane,sampleid,gerald_dir,demultiplexed_dir):
		"""Handles renaming a single samples sequence file or 
		files (if paired end) within a single bin."""
		# Identify the files to be processed.
		#pattern = "s_%d_([12]_)*sequence.txt" % lane
		#for entry in os.listdir(gerald_dir):
		#	if re.match(pattern,entry) is not None:
		#		filelist.append(entry)
		#filelist.sort()
		run=self.id

		# Generate list of files to be processed.
		if self.IsPaired:
			# Paired-end run
			filelist = [ "s_%d_1_sequence.txt" % lane, "s_%d_2_sequence.txt" % lane ]
		else:
			# Single-end run
			filelist = [ "s_%d_sequence.txt" % lane ]

		# Alter the sample id from nnnnLn to nnnnXn. Doing this so the file name
		# contains the same "X" notation of the samples in GNomEx.
		sampleid=re.sub("L","X",sampleid)

		retval = True
		for file in filelist:
			if self.IsPaired:
				end=file.split("_")[2]
			else:
				end=""
			orig_name=os.path.join(gerald_dir,file)
			if end:
				# paired-end
				new_name=os.path.join(demultiplexed_dir,"%s_%s_%s_%s.txt" % ( sampleid, run, lane, end ) )
			else:
				# Single end.
				new_name=os.path.join(demultiplexed_dir,"%s_%s_%s.txt" % ( sampleid, run, lane ) )
			self.datafiles.append(new_name)
			# Check if file already renamed and compressed.
			if os.path.exists(new_name+".gz"):
				self.Log(["File",orig_name,"already renamed and compressed."])
				continue
			if os.path.exists(orig_name):
				if os.path.exists(new_name):
					# Problem! Both the original or new file exist.
					self.Log(["PROBLEM! Found both",orig_name,"and",new_name])
					retval = False
				else:
					# Good! Rename the original file.
					self.Log(["Renaming",orig_name,"to",new_name])
					os.rename(orig_name,new_name)
			else:
				if os.path.exists(new_name):
					# OK - file has already been renamed.
					self.Log(["File",orig_name,"already renamed to",new_name])
				else:
					# Problem! File not found by either name.
					self.Log(["PROBLEM! Neither",orig_name,"nor",new_name,"found."])
					retval = False
		return retval

	def ProcessBin( self, bin, demultiplexed_dir ):
		"""Renames data files in a single bin subdirectory following demultiplexing.
		Data files will end up in demultiplexed_dir.
		"""
		os.chdir(bin)
		# Identify the GERALD directory.
		gerald_dir = None
		for entry in os.listdir("."):
			if os.path.isdir(entry) and entry[0:6] == 'GERALD':
				gerald_dir = entry
				break
		if gerald_dir is None:
			self.Log("PROBLEM! No GERALD directory in bin %s.\n" % bin )
			os.chdir('..')
			return False
		
		# Locate the SampleSheet.csv file.
		try:
			ifs = open("SampleSheet.csv")
			# Feedback.
			self.Log( "Processing SampleSheet.csv in %s, Gerald directory %s." % (bin, gerald_dir))
			# Read and discard header.
			rec = ifs.readline()
			# For each line...
			for rec in ifs:
				f = rec.strip().split(",")
				flowcell_barcode=f[0]
				lane=int(f[1])
				sampleid=f[2].upper()
				self.ProcessSample(lane,sampleid,gerald_dir,demultiplexed_dir)
			ifs.close()
		except IOError:
			# No sample sheet!
			self.Log("PROBLEM! No SampleSheet.csv file in bin %s.\n" % bin )
			os.chdir('..')
			return False
		os.chdir('..')
		return True

	def DemultiplexRenameDataFiles_18(self):
		"""Combines small gzipped fastq files produced by CASAVA 1.8
		into a single gzipped fastq file per sample / end. File names
		will be added to self.datafiles, for copying to repository.
		If sample belongs to Amplicon Express, however, the files are
		left as-is, and the entire directory of files is copied into
		the GNomEx repository.
		10/11/2013 - commenting out the Amplicon Express-specific
		code, because they are running some conventional samples
		at HCI. Brett Milash.
		"""
		# Read and parse the sample sheet. This gives lane, sample,
		# and request information.
		sample_info = []
		ifs = open(self.SampleSheetName())
		for rec in ifs:
			f = rec.strip().split(',')
			if f[0] == 'FCID':
				# Skip the header.
				continue
			lane = int(f[1])
			sample=f[2]
			request = f[9]
			requester = f[5]
			sample_info.append( (lane,sample,request, requester) )
		ifs.close()

		# Construct cat statements to build the result files.
		self.DetermineRunType()
		self.ClearJobs()
		for ( lane, sample, request, requester ) in sample_info:
			basename=os.path.join(self.dirname,"Unaligned")
			if AmpliconExpress.custom_handling and requester == AmpliconExpress.fullname:
				# Owner of all Amplicon Express samples.
				pass
				# Simply add the name of the result directory
				# for the lane into the self.datafiles list. 
				# This directory will get copied later.
				newname = os.path.join(basename,"Project_"+request)
				if newname not in self.datafiles:
					self.datafiles.append(newname)
			elif self.IsPaired:
				# Paired-end sequencing
				for end in [ 1, 2 ]:
					newname = "%s/%s_%s_%d_%d.txt.gz" % ( basename, sample, self.id, lane, end )
					self.datafiles.append(newname)
					if os.path.exists(newname):
						self.Log("%s already exists. Skipping." % newname)
						continue
					command = "cat %s/Project_%s/Sample_%s/%s_*_L00%d_R%d_*.fastq.gz > %s" % ( basename, request, sample, sample, lane, end, newname )
					self.AddJob( command )
					self.Log(command)
			else:
				# Single-end sequencing.
				end = 1
				newname = "%s/%s_%s_%d.txt.gz" % ( basename, sample, self.id, lane )
				self.datafiles.append(newname)
				if os.path.exists(newname):
					self.Log("%s already exists. Skipping." % newname)
					continue
				command = "cat %s/Project_%s/Sample_%s/%s_*_L00%d_R%d_*.fastq.gz > %s" % ( basename, request, sample, sample, lane, end, newname )
				self.AddJob( command )
				self.Log(command)
		retval = self.RunJobs(4,verbose=True)

		return retval

	def GenerateChecksums( self ):
		"""
		Creates a .md5 checksum file for each gzipped fastq file.
		"""
		md5files = []
		# Generate MD5 checksum files for each data file.
		self.Log(["Generating MD5 checksum files for run",self.id])
		for filename in self.datafiles:
			# Some of the file names may be directories. Copy
			# them using cp -r. This is the case for requests
			# from Amplicon Express.
			if os.path.isdir(filename):
				pass
			else:
				# Generate command.
				directory=os.path.dirname(filename)
				fastq_file=os.path.basename(filename)
				md5file = fastq_file + ".md5"
				cmd="(cd %s; /usr/bin/md5sum %s > %s)" % ( directory, fastq_file, md5file )
				self.AddJob( cmd )
				md5files.append(os.path.join( directory, md5file) )
		self.RunJobs(5,verbose=True)
		# Add the md5 checksum files to the list of data files
		# to be distributed.
		self.datafiles += md5files
		return True

	def PatchPcrDistribute( self ):
		"""
		Copies data from Patch PCR runs to GNomEx.
		"""
		self.CopyDataFiles()
		return False

	def DistributeDemultiplexed_18( self ):
		# CASAVA 1.8 creates many gzipped fastq files per sample.
		# Cat these files together into single file per sample/end.
		if self.DemultiplexRenameDataFiles_18():
			# Generate the MD5 checksums.
			if self.GenerateChecksums():
				# Copy data files.
				if self.CopyDataFiles():
					# Notify that run is complete.
					self.NotifyComplete('')
					return True
		return False

	def DistributeDemultiplexed( self ):
		"""Compress data files from multiplexed samples and copy to server for distribution. 
		This requires a connection to GNomEx for sample information."""
		self.Log(["Distributing data files from run",self.id])

		report = ''
		if self.DemultiplexRenameDataFiles():
			if self.CompressDataFiles():
				if self.CopyDataFiles():
					# Notify that the run is complete and data available.
					self.NotifyComplete(report)
					return True
		return False

	def Gerald( self ):
		"""Runs the GERALD scripts to produce sequence files.
		If geraldConfig.txt not found, create it.
		If GERALD directory not found or geraldConfig.txt is newer,
			run GERALD.pl.
		Run make in Gerald directory.
		"""
		# Create gerald config file.
		gerald_config_file = os.path.join(self.dirname,"geraldConfig.txt")
		if not os.path.exists( gerald_config_file ):
			if not self.CreateGeraldConfig(gerald_config_file):
				return False
		# If GERALD directory does not already exist or is older
		# than config file, run GERALD.pl to create it.
		gerald_dir = self.FindGeraldDirectory()
		retval = 0
		if gerald_dir is None or os.path.getctime(gerald_dir) < os.path.getctime(gerald_config_file):
			self.Log("Running GERALD.pl on run %s." % self.id)
			child_dir = os.path.join(self.dirname,"Data","Intensities","BaseCalls")
			child_stdout=open(os.path.join(child_dir,"gerald.out"),'w')
			child_stderr=open(os.path.join(child_dir,"gerald.err"),'w')
			child_args=[os.path.join(params.casava_dir,"GERALD.pl"),os.path.join(self.dirname,"geraldConfig.txt"),"--EXPT_DIR=.","-make"]
			p = subprocess.Popen(args=child_args,cwd=child_dir,stdout=child_stdout,stderr=child_stderr)
			retval = p.wait()
			child_stdout.close()
			child_stderr.close()
			# If Gerald ran successfully, get the gerald directory name.
			if retval == 0:
				gerald_dir = self.FindGeraldDirectory()
		else:
			self.Log("GERALD.pl already run in %s." % self.id)

		# If GERALD run successful, or had been run before...
		if retval == 0:
			# Run make in the GERALD... directory.
			self.Log(["Running make (to create sequence files) on run",self.id])
			child_dir = gerald_dir
			child_stdout=open(os.path.join(child_dir,"make.out"),'w')
			child_stderr=open(os.path.join(child_dir,"make.err"),'w')
			child_args=["/usr/bin/make","-j","8"]
			p = subprocess.Popen(args=child_args,cwd=child_dir,stdout=child_stdout,stderr=child_stderr)
			retval = p.wait()
			child_stdout.close()
			child_stderr.close()
		return retval == 0
	
	def RenameDataFiles( self, lanes=[1,2,3,4,5,6,7,8],dest_dir=None,gerald_parent_dir=None ):
		"""Renames the sequence data files produced by a non-barcoded run. Returns True if successful."""
		self.Log(["Renaming data files for run",self.id])
		self.DetermineRunType()
		# Get list of sample names and their lanes.
		query = """select sample.number samplenum, flowcellchannel.number lanenum
	from flowcellchannel
	join sequencelane on sequencelane.idflowcellchannel = flowcellchannel.idflowcellchannel
	join sample on sample.idsample = sequencelane.idsample
	where flowcellchannel.filename = '%s'""" % self.id
		g = GNomEx.GNomExConnection()
		connection = g.GnConnect(params.db_user,params.db_password)
		c = connection.cursor()
		c.execute(query)
		results = c.fetchall()
		# Build mapping from s_<lane>_[<end>_]sequence.txt files to <sample>_<run>_<lane>.txt files.
		gerald_dir = self.FindGeraldDirectory(gerald_parent_dir)
		self.Log(["RenameDataFiles: gerald_parent_dir = %s, gerald_dir = %s." % ( gerald_parent_dir, gerald_dir )])
		if dest_dir is None:
			dest_dir=gerald_dir
		filenames = []
		# List of full pathnames of all data files produced by this run.
		for rec in results:
			sample = rec[0]
			lane = rec[1]
			# Skip lane if it is not in one of the lanes to be processed. By default all lanes
			# are processed - only exception is for non-barcoded lanes on a flow cell with 
			# some barcoded samples.
			if lane not in lanes:
				continue
			if self.IsPaired:
				for end in [ 1, 2 ]:
					orig_name = os.path.join(gerald_dir,"s_%d_%d_sequence.txt" % ( lane, end ))
					new_name = os.path.join(dest_dir,"%s_%s_%d_%d.txt" % ( sample, self.id, lane, end ))
					filenames.append((orig_name,new_name))
					self.datafiles.append(new_name)
			else:
				orig_name = os.path.join(gerald_dir,"s_%d_sequence.txt" % lane)
				new_name = os.path.join(dest_dir,"%s_%s_%d.txt" % ( sample, self.id, lane ))
				filenames.append((orig_name,new_name))
				self.datafiles.append(new_name)
		# Rename the files.
		retval = True
		for ( orig_name, new_name ) in filenames:
			# Check if file already renamed and compressed.
			if os.path.exists(new_name+".gz"):
				self.Log(["File",orig_name,"already renamed and compressed."])
				continue
			if os.path.exists(orig_name):
				if os.path.exists(new_name):
					# Problem! Both the original or new file exist.
					self.Log(["PROBLEM! Found both",orig_name,"and",new_name])
					retval = False
				else:
					# Good! Rename the original file.
					self.Log(["Renaming",orig_name,"to",new_name])
					os.rename(orig_name,new_name)
			else:
				if os.path.exists(new_name):
					# OK - file has already been renamed.
					self.Log(["File",orig_name,"already renamed to",new_name])
				else:
					# Problem! File not found by either name.
					self.Log(["PROBLEM! Neither",orig_name,"nor",new_name,"found."])
					retval = False
		return retval
			
	def CompressDataFiles( self ):
		"""Compresses data files from a non-barcoded run. Returns True if successful."""
		# The list self.datafiles was initialized by RenameDataFiles.
		self.Log(["Compressing data files from run",self.id])
		#self.Log(["Data files:",self.datafiles])
		for filename in self.datafiles:
			if os.path.exists(filename+".gz"):
				self.Log(["File",filename,"already compressed."])
			else:
				self.Log(["Compressing",filename])
				self.AddJob("/bin/gzip -f " + filename)
		self.RunJobs(7,verbose=True)
		# Rename the list of data files with the .gz extension.
		for i in range(0,len(self.datafiles)):
			self.datafiles[i] += ".gz"
		return True

	def CopyDataFiles( self ):
		"""Copies data files from a run to the directories accessed by GNomEx."""
		# Get list of sample numbers, their request numbers, and the year of the request.
		query = """select distinct sample.number samplenum, request.number reqnum, request.createDate reqdate
		from flowcellchannel
		join sequencelane on sequencelane.idflowcellchannel = flowcellchannel.idflowcellchannel
		join sample on sample.idsample = sequencelane.idsample
		join request on sample.idrequest = request.idrequest
		where flowcellchannel.filename = '%s'""" % self.id

		g = GNomEx.GNomExConnection()
		connection = g.GnConnect(params.db_user,params.db_password,asdict=True)
		c = connection.cursor()
		c.execute(query)
		results = c.fetchall()

		# Set the umask so the files are world readable (and directories are
		# readable and executable).
		os.umask(002)

		# Build a mapping from sample number to result directory.
		result_directory = {}
		for rec in results:
			sample = rec['samplenum']
			reqnum = rec['reqnum']
			reqdate = rec['reqdate']
			# Edit the request number to remove any trailing digits.
			reqnum=reqnum[0:reqnum.index('R')+1]
			reqdate = str(reqdate.year)
			resultdir = os.path.join(params.repository_data_root[self.corefacilityname],reqdate,reqnum,"Fastq")
			result_directory[sample] = resultdir
			# Also index the result directory by request number.
			# This is used when copying directories.
			result_directory[rec['reqnum']] = resultdir

			# Create the result directory if necessary.
			if not os.path.exists( resultdir ):
				self.Log(["Making directory",resultdir])
				os.mkdir(resultdir)
			else:
				self.Log(["Result directory",resultdir,"already exists."])

		# Set up list of copy commands.
		for filename in self.datafiles:
			# Some of the file names may be directories. Copy
			# them using cp -r. This is the case for requests
			# from Amplicon Express.
			if os.path.isdir(filename):
				request_num = os.path.basename(filename).split('_')[-1]
				try:
					resultdir = result_directory[request_num]
					self.Log(["Copying",filename,"to",resultdir])
					#self.AddJob("/bin/cp -r " + filename + " " +resultdir)
					self.AddJob("/usr/bin/rsync " + filename + " " +resultdir)
				except KeyError:
					self.Log(["WARNING - request",request_num,"produced a data directory, but not found in database. Unable to copy file to result directory."])
			else:
				# Regular file, not a directory.
				sample_num = os.path.basename(filename).split('_')[0]
				try:
					resultdir = result_directory[sample_num]
					self.Log(["Copying",filename,"to",resultdir])
					# This command will loop until file is copied successfully. Brett Milash 8/6/2013.
					if filename.endswith(".gz"):
						self.AddJob( "until /bin/gunzip -c %s/%s > /dev/null; do /usr/bin/rsync %s %s; done" % ( resultdir, os.path.basename(filename), filename, resultdir ) )
					else:
						self.AddJob("/bin/cp " + filename + " " +resultdir)
				except KeyError:
					self.Log(["WARNING - sample",sample_num,"produced a data file, but not found in database. Unable to copy file to result directory."])
		self.RunJobs(6,verbose=True)

		return True

	def NotifyStarting( self ):
		"""Sends an email that the pipeline software is starting."""
		to_addr = params.lab_staff_addresses[self.corefacilityname]
		subject = "HiSeq pipeline starting. (%s)" % self.id
		message = "The data transfer for run %s is complete and the HiSeq pipeline is starting to process the run.\n" % self.id
		self.send_msg( self.FromAddr(), to_addr, subject, message )

	def NotifyComplete( self, run_summary_report='' ):
		"""Sends an email that the pipeline software has finished running, and that
		the data files are (or should be) on the server."""
		# Notify that the run is complete and data available.
		to_addr = params.lab_staff_addresses[self.corefacilityname]
		subject = "HiSeq pipeline processing is complete. (%s)" % self.id
		message = "The HiSeq pipeline has finished processing run %s.\n" % self.id
		message += "Please confirm that the sequence files are available for download\n"
		message += "before marking the run complete in GNomEx.\n"
		message += "The run folder will now be cleaned up to prepare it for archiving.\n"
		# Append run summary report, if any.
		if run_summary_report:
			message += '\n' + run_summary_report
		self.send_msg( self.FromAddr(), to_addr, subject, message )

	def Distribute( self ):
		"""Compress data files and copy to server for distribution. This requires
		a connection to GNomEx for sample information."""
		self.Log(["Distributing data files from run",self.id])
		if self.RenameDataFiles():
			if self.CompressDataFiles():
				if self.CopyDataFiles():
					# Notify that the run is complete and data available.
					self.NotifyComplete()
					return True
		return False
	
	def UseBasesMask( self, barcode_len, is_dual_index ):
		"""
		Returns the use_bases_mask value for the Illumina pipeline,
		which depends upon whether the flow cell has 1 or 2 data
		reads, whether the flow cell has 1 or 2 index reads, and whether
		the lane uses any of the index reads.
		"""
		# Every run has at least one data read.
		use_bases_mask = 'Y*'
		# Determine the mask given the barcode length.
		if barcode_len > 0:
			barcode_mask = ',I%dn*' % barcode_len
		else:
			barcode_mask = ',n*'

		# If there are any index reads, specify the first one.
		if self.num_index_reads > 0:
			use_bases_mask += barcode_mask
		# If there's a second index read...
		if self.num_index_reads == 2:
			# ... and this lane uses it...
			if is_dual_index:
				# ... add the mask for it.
				use_bases_mask += barcode_mask
			else:
				# ... else ignore it.
				use_bases_mask += ',n*'

		# If run has a second data read, use it.
		if self.num_data_reads == 2:
			use_bases_mask += ',Y*'

		return use_bases_mask


	def BclConvertComplex( self ):
		"""If flow cells contains lanes with varying lengths of barcode sequence, need to run
		pipeline seperately for each group of lanes with same lengths of barcode."""
		if self.sample_sheet is None:
			# Try to find the sample sheet. It should exist.
			self.sample_sheet = self.SampleSheetName()
		# Read the sample sheet and split it into separate sheets, each one of which has the lanes with 
		# barcodes of the same length.
		ifs = open(self.sample_sheet)
		lines = ifs.readlines()
		ifs.close()
		ofs = {}
		samplesheets = []

		for line in lines[1:]:
			# Get the length of the barcode for this line.
			f = line.split(',')
			barcode=f[4]
			barcode_len=len(barcode)
			# Look for a hyphen in the barcode. If found, this is
			# dual-indexed sample.
			is_dual_index=barcode.find('-')>-1

			# Write the line to the appropriate output file.
			try:
				ofs[barcode_len].write(line)
			except KeyError:
				# Generate file name.
				fname = re.sub(".csv","_%d.csv"%barcode_len,self.sample_sheet)
				samplesheets.append( (barcode_len,fname,is_dual_index) )
				ofs[barcode_len] = open(fname,'w')
				ofs[barcode_len].write(lines[0])	# Header.
				ofs[barcode_len].write(line)
		# Close the new sample sheet files.
		for filehandle in ofs.values():
			filehandle.close()

		# Process each sample sheet.
		retvals = []
		self.DetermineRunType()
		for (barcode_len, fname, is_dual_index) in samplesheets:
			# Determine the number of mismatches to allow.
			# For now this is based on barcode length.
			num_mismatches=1

			# No longer processing 8-base barcodes with 0 
			# mismatches. This was put in place for Scott Watkins
			# in the Jorde lab, but it is interferring with other
			# groups, and the Jorde lab just wants all the
			# reads anyway. Brett Milash, 9/26/2013.
			#if barcode_len == 8:
			#	num_mismatches=0
			
			# Make the full path of the sample sheet.
			sample_sheet_full_path=os.path.join(self.dirname,"Data","Intensities","BaseCalls",fname)
			outdir = os.path.join(self.dirname,"Unaligned_%d" % barcode_len)
			self.output_dirs.append( outdir )
			self.Log(["Starting bcl conversion, barcode length",barcode_len,", sample sheet", sample_sheet_full_path,", output dir", outdir])


			# If lane has dual index, recalculate the barcode_len
			# as the length of each of the two indices. Assuming 
			# they will be the same length as each other.
			if is_dual_index:
				barcode_len = (barcode_len-1)/2

			#if self.IsPaired:
			#	if is_dual_index:
			#		use_bases_mask = "Y*n,I%dn,I%dn,Y*n" % (barcode_len,barcode_len)
			#	elif barcode_len > 0:
			#		use_bases_mask = "Y*n,I%dn*,Y*n" % barcode_len
			#	else:
			#		use_bases_mask = "Y*n,n*,Y*n"
			#else:
			#	if is_dual_index:
			#		use_bases_mask = "Y*n,I%dn,I%dn" % (barcode_len,barcode_len)
			#	elif barcode_len > 0:
			#		use_bases_mask = "Y*n,I%dn" % barcode_len
			#	else:
			#		use_bases_mask = "Y*n,n*"
			use_bases_mask = self.UseBasesMask( barcode_len, is_dual_index )
			result = self.BclConvert( sample_sheet_full_path, outdir, use_bases_mask, compress_bcls=False, mismatches=num_mismatches )
			retvals.append( result )
			if not result:
				self.Log(["Bcl conversion of samples with barcode length", barcode_len, "failed."])
				break
			else:
				self.Log(["Bcl conversion of samples with barcode length", barcode_len, "succeeded."])
		# A side effect of calling BclConvert() is to replace
		# self.sample_sheet with one of the barcode-length-specific
		# sample sheets. This causes the problem when running the QC
		# report of reporting on only a subset of the lanes. By
		# restoring the self.sample_sheet we should get around the
		# QC report problem. Brett Milash 06/17/2014.
		self.sample_sheet = self.SampleSheetName()

		# If all Bcl file conversions successful, compress the
		# bcl files.
		if min(retvals) == True:
			self.CompressFiles("*.bcl")
		return min(retvals)

	def MergeBclConversions( self ):
		"""Combines results of multiple Bcl conversions (due to different length barcodes) into a single directory."""
		# Check if "Unaligned" directory exists. If not, create it.
		unaligned_dir = os.path.join(self.dirname, "Unaligned")
		if not os.path.exists( unaligned_dir ):
			os.mkdir( unaligned_dir )
		# Check if Undetermined_indices directory exists. If not, 
		# create it.
		undetermined_dir = os.path.join( unaligned_dir, "Undetermined_indices" )
		if not os.path.exists( undetermined_dir ):
			os.mkdir( undetermined_dir )

		# Process each output directory created during Bcl file conversion.
		if not self.output_dirs:
			# Find the output directories.
			for file in os.listdir(self.dirname):
				if file[0:10] == 'Unaligned_':
					self.output_dirs.append(file)
		for outputdir in self.output_dirs:
			outputdir = os.path.join(self.dirname,outputdir)
			# Get list of files in outputdir. Move any named "Project_*" 
			# to the unaligned directory.
			pattern = "Project_.*"
			for file in os.listdir( outputdir ):
				if re.match(pattern,file):
					try:
						os.rename( os.path.join( outputdir, file ), os.path.join( unaligned_dir, file ) )
					except OSError:
						# Problem! Folder already exists, probably from another lane.
						# mv Unaligned_%d/Project_?/Sample_* Unaligned/Project_?
						self.Log("Unaligned/%s directory already exists. Moving sample directories." % file)
						command="mv %s/%s/Sample_* %s/%s" % ( outputdir, file, unaligned_dir, file )
						self.Log(command)
						os.system(command)

						#message = "Problem! OSError generated trying to move %s from %s to %s." % ( file,outputdir,unaligned_dir)
						#self.Log( message )
						#self.Notify("Problem with flow cell processing for " + self.dirname, message )

			# Look for directories containing bad barcode reads. 
			# These are used by the QC report. If any found, move 
			# to the undetermined indices directory.
			bad_barcodes_dir = os.path.join( outputdir, "Undetermined_indices" )
			if os.path.exists( bad_barcodes_dir ):
				pattern = "Sample_lane.*"
				for file in os.listdir( bad_barcodes_dir ):
					if re.match(pattern,file):
						os.rename( os.path.join( bad_barcodes_dir, file ), os.path.join( undetermined_dir, file ) )
				
		return True

	def Qc_18( self ):
		"""Runs a barcode QC report on a flowcell for pipeline version
		1.8."""
		return self.Qc( pipeline_version="1.8" )

	def FindQcLanes( self ):
		"""
		Finds lanes in the current run that should show up on the QC
		report.
		"""
		# Open the sample sheet.
		if self.sample_sheet is None:
			# Try to find the sample sheet. It should exist.
			self.sample_sheet = self.SampleSheetName()
		lanes = set()
		ifs = open(self.sample_sheet)
		for rec in ifs:
			f = rec.strip().split(',')
			if f[0] == 'FCID':
				# Skip the header.
				continue
			lane = int(f[1])
			requester = f[5]
			if AmpliconExpress.custom_handling and requester == AmpliconExpress.fullname:
				continue
			else:
				lanes.add(lane)
		ifs.close()
		l = list(lanes)
		l.sort()
		return l

	def FindQcSkippedLanes( self ):
		"""Identify any lanes that we should skip in the QC report."""
		# Open the sample sheet.
		if self.sample_sheet is None:
			# Try to find the sample sheet. It should exist.
			self.sample_sheet = self.SampleSheetName()
		skip_lanes = set()
		ifs = open(self.sample_sheet)
		for rec in ifs:
			f = rec.strip().split(',')
			if f[0] == 'FCID':
				# Skip the header.
				continue
			lane = int(f[1])
			requester = f[5]
			if AmpliconExpress.custom_handling and requester == AmpliconExpress.fullname:
				skip_lanes.add(lane)
		ifs.close()
		return skip_lanes

	def Qc( self, pipeline_version="1.7" ):
		"""Runs a barcode QC report on the flowcell and sends the
		results via email. Also saves a copy to the FlowCellData directory
		for this flow cell."""
		try:
			# Identify any lanes that we should skip.
			#lanes=range(1,9)
			#skip_lanes = self.FindQcSkippedLanes()
			#for l in skip_lanes:
			#	lanes.remove(l)

			# Identify the lanes that should show up on the QC
			# report. This approach handles regular flow cells
			# and rapid-run 2 lane flow cells. It also handles
			# Amplicon Express requests. Brett Milash, 10/11/2013.
			lanes = self.FindQcLanes()

			# Generate the report.
			import EvaluateIndexing
			outputfile=os.path.join(self.dirname,"barcode_report_%s.xls"% self.id )
			self.Log(["About to run qc report for",self.id,"lanes",`lanes`,"writing to",outputfile])
			EvaluateIndexing.RunReport( self.dirname, outputfile, pipeline_version, lanes )

			# Locate the flowcelldata directory for this flow cell.
			# The flowcell barcode is the last part of the "_" delimited
			# self.id, following the first letter (which is the flow cell 
			# position in the sequencer, A or B.
			barcode = self.id.split("_")[-1]
			# Trim the initial character from the run folder barcode if this is a 
			# hiseq 2000 run.
			if barcode[0] in ['A','B']:
				barcode = barcode[1:]
			g = GNomEx.GNomExConnection()
			connection = g.GnConnect(params.db_user,params.db_password,asdict=True)
			c = connection.cursor()
			c.execute("select number, createdate from flowcell where barcode = %s", barcode)
			results = c.fetchall()
			connection.close()
			fcnumber=results[0]["number"]
			fcdate=results[0]["createdate"]
			fcyear = str(fcdate.year)
			flowcelldatadir=os.path.join(params.repository_root_dir,"FlowCellData",fcyear,fcnumber)

			# Save the report to that directory.
			self.Log(["Saving qc report for",self.id,"to",flowcelldatadir])
			os.system("cp "+outputfile+" "+flowcelldatadir)

			# Email the report.
			to_addr = params.lab_staff_addresses[self.corefacilityname]
			subject = "HiSeq pipeline QC report (%s)" % self.id
			message = "Here is the barcode processing report for run %s." % self.id
			self.Log(["Sending qc report for",self.id,"to"]+to_addr)
			self.send_msg( self.FromAddr(), to_addr, subject, message, attachments=[outputfile] )
			return True
		except:
			self.Log(["PROBLEM! Can't generate/send QC report."])
			raise

	def RunShellCommand( self, command, pattern, numjobs=8 ):
		"""Runs command on all files matching pattern. Runs in 
		parallel numjobs at a time, defaults to 8 concurrent jobs."""

		self.Log(["Running command", command, "on files matching pattern",pattern,"from run",self.id])
		cmd = '/usr/bin/find %s -name "%s" -print' % ( self.dirname, pattern )
		child_dir = self.dirname
		child_stderr=open(os.path.join(child_dir,"find.err"),'a')
		p = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE).stdout
		for rec in p:
			fname = rec.strip()
			cmd = "%s '%s'" % (command,fname)
			self.AddJob(cmd)
		p.close()
		self.RunJobs(numjobs,verbose=True)

	def CompressFiles( self, pattern, numjobs=8 ):
		"""Compresses files under run folder that match the given pattern."""
		self.Log(["Compressing files matching",pattern,"from run",self.id])
		self.RunShellCommand("/bin/gzip -f", pattern, numjobs )

	def UncompressFiles( self, pattern, numjobs=8 ):
		"""Unompresses files under run folder that match the given pattern."""
		self.Log(["Unompressing files matching",pattern,"from run",self.id])
		self.RunShellCommand("/bin/gunzip", pattern, numjobs )
	

	def Cleanup( self ):
		"""Clean up run folder in preparation for archiving. Qseq files
		removed, text files and bcl files compressed."""
		self.Log(["Cleaning up run",self.id])
		# Remove qseq files.
		self.Log(["Removing qseq files from run",self.id])
		child_dir = self.dirname
		child_stdout=open(os.path.join(child_dir,"find.out"),'a')
		child_stderr=open(os.path.join(child_dir,"find.err"),'a')
		child_args=["find",self.dirname,"-name","*qseq.txt","-delete"]
		p = subprocess.Popen(args=child_args,cwd=child_dir,stdout=child_stdout,stderr=child_stderr)
		retval = p.wait()
		child_stdout.close()
		child_stderr.close()

		# Compress .bcl files.
		# Now doing this right after BCL to QSEQ conversion.
		#self.CompressFiles("*.bcl")

		# Compress .txt files.
		self.CompressFiles("*.txt")

		# Compress .control files.
		self.CompressFiles("*.control")
		return retval == 0

	def Archive( self ):
		"""Run QC report and notify someone that this run is ready for archiving."""
		to_addr = params.archive_addresses
		subject = "HiSeq Pipeline Notification"
		message = "Run %s is complete, cleaned up, and ready for archiving." % self.id
		self.Log(["Notifying",to_addr,"that",message])
		try:
			self.send_msg( self.FromAddr(), to_addr, subject, message, attachments=[] )
			return True
		except:
			self.Log(["PROBLEM! Can't notify",to_addr,"that",message])
			return False
	
	def SampleSheetName( self ):
		"""Returns full-path name of samplesheet file."""
		return os.path.join(self.dirname,"Data","Intensities","BaseCalls","created_samplesheet.csv")

	def Reprocess( self ):
		"""Uncompress bcl and txt files so run can be processed again from the beginning."""
		self.Log(["Reprocessing run", self.id, "run directory", self.dirname])
		# Remove the sample sheet file, to force one to be regenerated.
		#samplesheet_fname = self.SampleSheetName()
		#if os.path.exists(samplesheet_fname):
		#	os.unlink(samplesheet_fname)

		# Rename the unaligned directory, if any, so a new one
		# can be created. Pipeline won't overwrite this directory
		# unless told to do so.
		unaligned_dir=os.path.join(self.dirname,"Unaligned")
		if os.path.exists( unaligned_dir ):
			i = 1
			newname = '%s.%d' % (unaligned_dir,i)
			# Increment i until we find a name that doesn't
			# exist yet.
			while os.path.exists( newname ):
				i+=1
				newname = '%s.%d' % (unaligned_dir,i)
			self.Log(["Renaming",unaligned_dir,"to",newname])
			os.rename(unaligned_dir,newname)

		# Uncompress the bcl files.
		flowcell = self.FlowCellVersion()
		if flowcell == "HiSeq Flow Cell v3":
			self.UncompressFiles("*.bcl.gz")
		# Uncompress the txt files.
		self.UncompressFiles("*.txt.gz")
		# Uncompress the control files.
		self.UncompressFiles("*.control.gz")
		return True
	
	def CheckIfMiseq( self ):
		"""
		Checks if run is from a Miseq sequencer.
		"""
		self.Log("Checking if run %s is a miseq run." % self.id)
		configfile = os.path.join(self.dirname, "RunInfo.xml")
		doc = xml.dom.minidom.parse(configfile)
		e=doc.getElementsByTagName("Instrument")[0]
		instrument=e.firstChild.nodeValue
		is_miseq = instrument in params.miseq_ids
		if is_miseq:
			self.Log("Run %s is a miseq run." % self.id )
		else:
			self.Log("Run %s is a hiseq run." % self.id )
		return is_miseq
	
	def CheckIfPatchPcr( self ):
		"""
		Determines if a miseq run is a patch pcr run.
		"""
		self.Log("Checking if run %s is a patch pcr run." % self.id)
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
		else:
			self.Log("Run %s is not a patch pcr run." % self.id)
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
			"--sample-sheet",self.SampleSheetName(),
			"--barcode-mismatches","1",
			"--use-bases-mask",use_bases_mask,
			"--minimum-trimmed-read-length", `min_read_length`
		]
		self.Log("PatchPcrDemultiplex executing command: %s" % `child_args`)

		p = subprocess.Popen(args=child_args,cwd=child_dir,stdout=child_stdout,stderr=child_stderr)
		retval = p.wait()
		child_stdout.close()
		child_stderr.close()
		if retval == 0:
			return True
		else:
			message="Demultiplexing of patch PCR run failed. See bcl2fastq.err in %s for details." % self.dirname
			self.Log( message )
			self.Notify( "Run %s demultiplexing failed" % self.dirname, message )
			return False
	
	def GetSampleSheetProjectsSamples(self):
		"""
		Reads the sample sheet for the run and returns a dictionary
		indexed by project of all the sample ids.
		"""
		d = {}
		# Open the sample sheet.
		ifs=open(self.SampleSheetName())
		# Read until the header is found.
		for rec in ifs:
			if rec.startswith("Sample_ID"):
				break

		samplenum=0
		for rec in ifs:
			samplenum+=1
			f=rec.strip().split(',')
			project=f[8]
			# Illumina adds _S# to the sample id based on its
			# position in the sample sheet.
			sample=f[0]+"_S%d"%samplenum
			try:
				d[project].append(sample)
			except KeyError:
				d[project] = [sample]
		ifs.close()
		self.Log("GetSampleSheetProjectsSamples: returning %s." % `d` )
		return d

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
		r1_ifs=Fastq.Reader(r1_in_file)
		r2_ifs=Fastq.Reader(r2_in_file)
		r1_ofs=Fastq.Writer(r1_out_file)
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
			r3_ifs=Fastq.Reader(r3_in_file)
			r2_ifs=Fastq.Reader(r2_in_file)
			r2_ofs=Fastq.Writer(r2_out_file)
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
				self.Log("Postprocessing project %s sample %s." % ( project,sample))
				self.PatchPcrPostProcessSample(project,sample)
	
		# Generate the MD5 checksums.
		return self.GenerateChecksums()

	def DbAdd( self, dbconnection ):
		"""Insert object into database."""
		c = dbconnection.cursor()
		try:
			c.execute("insert into run (id,dirname,state) values (%s,%s,%s)", (self.id,self.dirname,self.state))
		except TypeError:
			print 'self.id',self.id
			print 'self.dirname',self.dirname
			print 'self.state',self.state
			raise
		dbconnection.commit()

	def DbUpdate( self, dbconnection ):
		c = dbconnection.cursor()
		c.execute("update run set state=%s where id=%s", (self.state,self.id,))
		dbconnection.commit()

class PipelineMgr(Logger,Emailer):
	"""PipelineMgr manages running the HiSeq Pipeline."""
	def __init__( self, runs_directories ):
		self.db_file = "pipelinemgr.db"

		# List of parent directories in which runs folders are created.
		self.runs_directories = runs_directories

		# Verify that runs directories exist.
		for runs_dir in self.runs_directories:
			if not os.path.exists(runs_dir):
				raise IOError, "Runs directory '%s does not exist." % runs_dir

		# Open PipelineMgr database that tracks processing of runs. If it doesn't
		# exist yet, create it.
		if self.DbExists():
			self.DbOpen()
		else:
			self.DbCreate()

		# Define state transitions.
		self.transition = {
			# curr_state, function, success_state, fail_state
			States.new: (Run.CheckRegisteredVerbose,
				States.wait_for_data,
				States.check_registered ),
			States.check_registered: (Run.CheckRegisteredSilent,
				States.wait_for_data,
				States.check_registered),
			States.wait_for_data: (Run.CheckTransferComplete,
				States.check_if_miseq,
				States.wait_for_data),
			States.make_sample_sheet: (Run.MakeSampleSheet,
				States.check_single_index_length,
				States.error_detected),
			States.check_single_index_length: (Run.CheckSingleIndexLength,
				States.simple_bcl_convert,
				States.complex_bcl_convert),
			States.simple_bcl_convert: (Run.BclConvert,
				States.distribute_demultiplexed,
				States.error_detected),
			States.complex_bcl_convert: (Run.BclConvertComplex,
				States.merge_bcl_conversions,
				States.error_detected),
			States.merge_bcl_conversions: (Run.MergeBclConversions,
				States.distribute_demultiplexed,
				States.error_detected),
			States.distribute_demultiplexed: (Run.DistributeDemultiplexed_18,
				States.cleanup,
				States.error_detected),
			States.cleanup: (Run.Cleanup,
				States.qc,
				States.error_detected),
			States.qc: (Run.Qc_18,
				States.archive,
				States.error_detected),
			States.archive: (Run.Archive,
				States.complete,
				States.error_detected),
			States.error_detected: (Run.NotifyError,
				States.error_notified,
				States.error_detected),
			States.reprocess: (Run.Reprocess,
				States.make_sample_sheet,
				States.error_detected ),
			States.check_if_miseq: (Run.CheckIfMiseq,
				States.check_if_patchpcr,
				States.make_sample_sheet),
			States.check_if_patchpcr: (Run.CheckIfPatchPcr,
				States.patchpcr_sample_sheet,
				States.complete),
			States.patchpcr_sample_sheet: (Run.PatchPcrSampleSheet,
				States.patchpcr_demultiplex,
				States.error_detected),
			States.patchpcr_demultiplex: (Run.PatchPcrDemultiplex,
				States.patchpcr_postprocess,
				States.error_detected),
			States.patchpcr_postprocess: (Run.PatchPcrPostprocess,
				States.patchpcr_distribute,
				States.error_detected),
			States.patchpcr_distribute: (Run.PatchPcrDistribute,
				States.qc,
				States.error_detected),
		}

	def DbExists( self ):
		"""Checks if PipelineMgr database file exists."""
		return os.path.exists(self.db_file)
	
	def DbCreate( self ):
		"Create and open the database."""
		self.Log("Creating database.")
		self.DbOpen()
		c = self.db_connection.cursor()

		# Create runs table.
		c.execute("create table run (id primary key,dirname,state)")
		self.db_connection.commit()

		c.close()
		self.Log("Database creation complete.")

	def DbOpen( self ):
		"Open the database."""
		#self.Log("Opening database.")
		self.db_connection = sqlite3.connect(self.db_file)
		
	def CleanUp( self ):
		"""Clean any old and deleted runs out of the database. Any
		Run in the database that doesn't appear on the filesystem
		can be removed from the database."""
		# Get all the runs out of the database.
		runs = self.GetAllRuns()
		c = self.db_connection.cursor()
		# For each run...
		for (id,dirname,state) in runs:
			# ... if not present on the file system...
			if not os.path.exists(dirname):
				# ... remove the record from the database.
				self.Log("Deleting old run %s (%s) from database." % (id,dirname))
				c.execute("delete from run where id = %s",id)
				self.db_connection.commit()

	def ProcessRuns( self ):
		# Discover any new runs.
		self.Log("Discovering new runs.")
		self.Discover()
		self.Log("Discovery complete.")

		# Select any runs from the database that are not in the 
		# done or error states, and sort them by (decreasing state,
		# increasing create date. This will prioritize so that the oldest
		# and furthest along will get priority.
		runs = self.GetActiveRuns()
		runs.sort()
		# For each run in the list...
		for run in runs:
			prev_state = None
			while run.state != prev_state:
				prev_state = run.state
				try:
					if self.transition.has_key(run.state):
						(method,success_state,fail_state) = self.transition[run.state]
						if method(run):
							# Function successful. 
							# Call state-specific methods for each lane.
							# Move run to new state.
							self.Log(["Run",run.id,"transitioned from",run.state,"to",success_state])
							run.state = success_state
						else:
							self.Log(["Run",run.id,"transitioned from",run.state,"to",fail_state])
							run.state = fail_state
						run.DbUpdate(self.db_connection)
					else:
						self.Log(["Skipping run",run.id,", no transitions from state", run.state])
				except:
					(etype,evalue,etraceback) = sys.exc_info()
					details = traceback.format_exception( etype,evalue,etraceback)
					self.Log(["Exception encountered processing",run.id,". Continuing with next run. Details:" ]+details)
				
		# Clean any old runs out of the database.
		self.CleanUp()

	def GetActiveRuns( self ):
		"""Returns list of active runs from database."""
		c = self.db_connection.cursor()
		c.execute("select id, dirname, state from run where state != %s and state != %s",(States.error_notified,States.complete))
		recs = c.fetchall()
		runs = []
		for rec in recs:
			runs.append(Run(rec[0],rec[1],rec[2]))
		return runs
	
	def GetAllRuns( self ):
		"""Return list of all known runs."""
		c = self.db_connection.cursor()
		c.execute("select id,dirname,state from run")
		recs = c.fetchall()
		return recs
	
	def GetAllRunIds( self ):
		"""Return list of all known run ids from database."""
		c = self.db_connection.cursor()
		c.execute("select id from run")
		recs = c.fetchall()
		runs = map( lambda x:x[0], recs )
		return runs

	def Discover( self ):
		"""Look for any new runs. Enter them into the database
		in the starting state."""
		# Find subdirectories in the runs directory whose name matches the
		# regular expression.
		known_run_ids = self.GetAllRunIds()
		#self.Log(["Known run ids",known_run_ids])

		new_runs = []
		# Pattern consists of Date_SequencerName_RunNumber_Flowcell.
		# Pattern updated to match both HiSeq and MiSeq runs.
		# Brett Milash 3/20/2015.

		# Use this pattern to find both Hiseq and Miseq runs:
		pattern = re.compile("[0-9]*_[A-Z0-9]*_[0-9]*_[A-Z0-9-]*$")

		# Use this pattern to find only Hiseq runs:
		#pattern = re.compile("[0-9]*_[A-Z0-9]*_[0-9]*_[AB][A-Z0-9]*")

		for runs_dir in self.runs_directories:
			for subdir in os.listdir(runs_dir):
				run_full_path=os.path.join(runs_dir,subdir)
				if not os.path.isdir(run_full_path):
					# Skip non-directory files.
					continue
				if pattern.match(subdir):
					#self.Log(["Checking subdirectory",subdir])
					# Found a runs folder.
					if subdir not in known_run_ids:
						# Discovered a new run folder. 
						new_runs.append((subdir,run_full_path))
					else:
						self.Log(["Subdirectory",subdir,"already known."])
		if new_runs:
			for (subdir,run_full_path) in new_runs:
				run = Run( subdir,run_full_path )
				run.DbAdd( self.db_connection )
				self.Log(["Discovered run",subdir])
		else:
			self.Log("No new sequencing runs.")

def main():
	# Set up semaphore to guarantee exclusivity.
	keyval = 42
	flags = sysv_ipc.IPC_CREX
	initial = 1
	try:
		sem = sysv_ipc.Semaphore(key=keyval, flags=flags, initial_value=initial)
		sem.acquire(timeout=0)
	except sysv_ipc.ExistentialError:
		Logger().Log( "%s - program already running." % sys.argv[0] )
		sys.exit(0)

	try:
		mgr = PipelineMgr( params.root_directories )
		mgr.ProcessRuns()
	except IOError, message:
		Logger().Log("%s - a runs directory in %s doesn't exist." % ( sys.argv[0], params.root_directories ))
	except:
		# Clean up after unforseen errors.
		Logger().Log("%s - unexpected exception caught." % sys.argv[0])
		sem.release()
		sem.remove()
		raise

	# Release and remove the semaphore.
	sem.release()
	sem.remove()

if __name__ == "__main__":
	main()
