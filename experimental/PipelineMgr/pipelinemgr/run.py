import os
import sys
import traceback
import socket
import xml.dom.minidom
import re

from . import db
from .state import State
from .emailer import Emailer
from .logger import Logger
from .simultaneousjobrunner import SimultaneousJobRunner
from . import pipelineparams
from .gnomex import GNomExConnection
from .util import EraseCommas, WriteSamplesheetHeader
from .runprocessor import CreateRunProcessor

class Run(Emailer,SimultaneousJobRunner):
	"""
	Run object represents one run folder. 
	Run objects store their ids, run folders, and processing state in
	a SQLite database to maintain state information when the 
	pipeline isn't running.
	"""
	def __init__(self,id,run_directory,state=State.new):
		self.id=id
		self.run_directory=run_directory
		self.state=state
		self.corefacilityname = None
		self.gnomex_connection = None
		self.GetRunDetails()
		self.done_file=os.path.join(self.run_directory,"Basecalling_Netcopy_complete_Read%d.txt" % self.run_num_reads)
		self.gnomex_connection = GNomExConnection().GnConnect(pipelineparams.db_user,pipelineparams.db_password)
		self.InitializeCoreFacility()
	
	def GetRunDetails(self):
		"""
		Determines whether sequencing run was from a hiseq or miseq sequencer,
		how many data reads (single-end vs paired-end) and how many index
		reads (0, 1, or 2) the run had. 
		Initializes instrument, data_reads, index_reads, and run_data_reads
		values.
		"""
		doc = xml.dom.minidom.parse(os.path.join(self.run_directory,"RunInfo.xml"))
		layout=doc.getElementsByTagName("FlowcellLayout")
		# Determine what instrument. Doing this by lane count, which
		# won't work if HiSeq instrument in rapid mode.
		if int(layout[0].getAttribute("LaneCount")) == 8:
			self.instrument='hiseq'
		else:
			self.instrument='miseq'

		# Determine the number of data reads and index reads.
		self.run_data_reads=0
		self.run_index_reads=0
		for tag in doc.getElementsByTagName("Read"):
			if tag.getAttribute("IsIndexedRead") == 'N':
				self.run_data_reads += 1
				self.run_data_read_length = int(tag.getAttribute("NumCycles"))
			elif tag.getAttribute("IsIndexedRead") == 'Y':
				self.run_index_reads += 1
				self.run_index_read_length = int(tag.getAttribute("NumCycles"))
		self.run_num_reads=self.run_data_reads+self.run_index_reads

	def FromAddr( self ):
		"""
		Returns from address for emailing purposes.
		"""
		# Unix account under which PipelineMgr.py runs.
		username = os.environ['LOGNAME']
		return username + "@" + socket.getfqdn()

	def NotifyLabStaff( self, subject, message ):
		"""
		Notify lab staff about a sequencing run issue.
		"""
		to_addr = pipelineparams.lab_staff_addresses
		self.send_msg( self.FromAddr(), to_addr, subject, message )
	
	def SmCheckRegisteredVerbose(self):
		if self.SmCheckRegisteredSilent():
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
	
	def Query( self, query ):
		"""
		Runs a SQL query against GNomEx, and retuns a cursor from
		which the results can be retrieved.
		"""
		c = self.gnomex_connection.cursor()
		c.execute(query)
		return c

	def SmCheckRegisteredSilent(self):
		"""
		This method determines if the run has been registered in
		GNomEx.
		"""
		#g = gnomex.GNomExConnection()
		#connection = g.GnConnect(pipelineparams.db_user,pipelineparams.db_password)
		#c = connection.cursor()
		#c.execute("select count(*) from flowcellchannel where filename = %s", self.id)
		c=self.Query("select count(*) from flowcellchannel where filename = '%s'"%self.id)
		
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
			if self.SmCheckTransferComplete():
				subject = "HiSeq Pipeline Problem! (%s)" % self.id
				message = "The sequencing run %s is done, but the run folder has not been entered in GNomEx.\n" % self.id
				message += "The pipeline can not start until the flow cell has been entered into GNomEx.\n"
				self.NotifyLabStaff( subject, message )
			return False
	
	def SmCheckTransferComplete(self):
		"""
		This method checks if the transfer of all files from the 
		sequencer is done.
		"""
		return os.path.exists(self.done_file) or os.path.exists(self.done_file+".gz")
	
	def SmRunPipeline(self):
		"""
		SmRunPipeline creates one or more RunProcessors to perform
		base calling, demultiplexing, and data distribution.
		"""
		# Create a sample sheet.
		samplesheet_fname = self.CreateSampleSheet()

		# Create RunProcessor(s) for run.
		runprocessors={}
		samplesheets={}
		sample_sheet_num=0

		# Read and parse the sample sheet file - doing it this way so
		# we can override what's in GNomEx by editing the file, rather
		# than relying on what's in GNomEx being correct. This works
		# well when reprocessing the pipeline.
		ifs=open(samplesheet_fname)

		# Read the discard the header.
		ifs.readline()

		# For each lane on flow cell get length of the index reads 
		# and the demultiplex method.
		for rec in ifs:
			(fcid,lane,sampleid,sampleref,index,description,control,demultiplex_method,Operator,SampleProject)=rec.strip().split(',')
			barcodes=index.split('-')
			try:
				barcode_a_len=len(re.sub("N","",barcodes[0]))
			except IndexError:
				barcode_a_len=0
			try:
				barcode_b_len=len(re.sub("N","",barcodes[1]))
			except IndexError:
				barcode_b_len=0

			key=(self.instrument,		# from RunInfo.xml
				self.run_data_reads,	# from RunInfo.xml
				self.run_index_reads,	# from RunInfo.xml
				barcode_a_len,		# from sample sheet
				barcode_b_len,		# from sample sheet
				demultiplex_method)	# from sample sheet
			# Create a new run processor. Also, create a new
			# Sample sheet for that run processor.
			if key not in runprocessors:
				self.Log("Need a new run processor for key %s." % str(key))
				# Open the sample sheet.
				sample_sheet_num+=1
				sheetname=os.path.join(self.run_directory,"created_samplesheet_%d.csv" % sample_sheet_num)
				self.Log("Creating sample sheet %s." % sheetname)
				samplesheets[key]=open(sheetname,'w')
				WriteSamplesheetHeader(samplesheets[key])
				# Create the run processor.
				rp=CreateRunProcessor(sample_sheet_num,self.run_directory,
					self.instrument,
					self.run_data_reads,
					self.run_index_reads,
					barcode_a_len,
					barcode_b_len,
					demultiplex_method,
					sheetname)
				runprocessors[key]=rp
			# Add the lane to the appropriate run processor
			runprocessors[key].AddLane(int(lane))
			samplesheets[key].write(rec)
		ifs.close()
		# Close all the new sample sheets.
		for ofs in samplesheets.values():
			ofs.close()

		self.Log("Created %d run processors for %s." % ( len(runprocessors),self.id))
		# Execute each RunProcessors' ProcessRun method.
		results=[]
		for runprocessor in list(runprocessors.values()):
			results.append(runprocessor.ProcessRun())
		self.Log("All run processors finished. Return values: %s." % str(results))

		# Collect the results, distribute to GNomEx.

		# Generate a QC report.

		# Clean up.
		#self.CompressFiles("*.bcl")
		#self.CompressFiles("*.txt")
		#self.CompressFiles("*.control")
		return True
	
	def SmReprocess(self):
		"""
		SmReprocess prepares an already-processed run for reprocessing.
		This is necessary if samples were mis-labeled in GNomEx for 
		example.
		"""
		# Generate a new sample sheet file with a different name. 
		# Identify lanes that have changed, and put those changed
		# lanes into the master sample sheet.

		# Uncompress any data files compressed after the initial 
		# processing.
		self.UncompressFiles("*.bcl")
		self.UncompressFiles("*.txt")
		self.UncompressFiles("*.control")
		return True

	def RunShellCommand( self, command, pattern, numjobs=8 ):
		"""
		Runs command on all files matching pattern. Runs in 
		parallel numjobs at a time, defaults to 8 concurrent jobs.
		Limited to files to Data subdirectory of run directory.
		"""
		self.Log(["Running command", command, "on files matching pattern",pattern,"from run",self.id])
		cmd = '/usr/bin/find %s -name "%s" -print' % ( os.path.join(self.run_directory,"Data"), pattern )
		child_dir = self.run_directory
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
		self.Log(["Uncompressing files matching",pattern,"from run",self.id])
		self.RunShellCommand("/bin/gunzip", pattern, numjobs )


	def DbAdd( self, dbconnection ):
		"""
		Inserts this Run object into database.
		"""
		c = dbconnection.cursor()
		try:
			c.execute("insert into run (id,run_directory,state) values (?,?,?)", (self.id,self.run_directory,self.state.value))
		except TypeError:
			print('self.id',self.id)
			print('self.run_directory',self.run_directory)
			print('self.state',self.state.value)
			raise
		dbconnection.commit()

	def DbUpdate( self, dbconnection ):
		"""
		Stores run's state in database.
		"""
		c = dbconnection.cursor()
		c.execute("update run set state=? where id=?", (self.state.value,self.id,))
		dbconnection.commit()
	
	def Process(self,dbconnection):
		"""
		This method uses the state machine to call various Run methods
		depending on the Run's current state.
		"""
		prev_state = None
		# While the run is still making progress from state to state ...
		while self.state != prev_state:
			prev_state = self.state
			try:
				# ... get the method to call for the run's current
				# state.
				(method,success_state,fail_state) = statemachine[self.state.value]
				# Set up a dictionary to map the method's return
				# value to the run's next state.
				next_state={
					True:success_state,
					False:fail_state,
				}
			
				try:
					# Call the method, and set the run's
					# new state depending on the method's
					# return value.
					self.state = next_state[ method(self) ]
					if self.state == prev_state:
						self.Log(["Finished with run",self.id,", no progress out of state",self.state.value])
					else:
						self.Log(["Run",self.id,"transitioned from",prev_state.value,"to",self.state.value])
					# Store the runs' new state.
					self.DbUpdate(dbconnection)
				except:
					# Handle exceptions thrown in method.
					(etype,evalue,etraceback) = sys.exc_info()
					details = traceback.format_exception( etype,evalue,etraceback)
					self.Log(["Exception encountered processing",self.id,". Continuing with next run. Details:" ]+details)
					break
			except KeyError:
				self.Log(["Finished with run",self.id,", no transitions from state", self.state.value])
				break
			except AttributeError:
				self.Log("Current state is '%s'." % self.state)
				raise

	def CreateSampleSheet(self,override_filename=None,lanes=[1,2,3,4,5,6,7,8]):
		"""
		Creates sample sheet for multiplexed flow cell, and writes 
		it to the BaseCalls directory within the run directory. 
		Returns full path of sample sheet file name.
		Changes for version 1.8 of CASAVA:
			Added SampleProject column
			altered other column names
			write file to BaseCalls directory
		Permits samples without barcodes. These samples will
			have an Null values for the sample.barcode field, and 
			no barcode processing will happen for that lane.
		If only single sample is in a lane, the barcode for
			that sample will be ignored, and it will be treated as 
			if a single non-barcoded sample is there. Doing this 
			because sequence quality for single barcodes is poor, 
			and many clients just want all the reads anyway.
		No longer supporting rapid runs. Just commented out the code.
			These runs are registered in GnomEx as a single lane, 
			but actually have 2 lanes of data (numbered 1 and 2). 
		The optional override_filename argument is used in flowcell
			reprocessing.

		Generate sample sheet name. If it already exists just
		return its name. Doing this to make it easier to override
		the automatic sample sheet creation. Brett Milash, 10/11/2013.
		Now putting these at the top level of the run directory,
		much easier to locate the files there. Brett Milash 06/01/2017.

		2017-07-10 - added pipeline protocol to sample sheet. If 
		pipeline protocol is undefined then 'None' is written to
		sample sheet. Putting this in the Recipe column, since that is 
		unused. Brett Milash.
		"""
		self.Log("Creating sample sheet for %s." % self.id)
		if override_filename is not None:
			samplesheet_fname = override_filename
		else:
			samplesheet_fname = os.path.join(self.run_directory,"created_samplesheet.csv")
		if os.path.exists(samplesheet_fname):
			return samplesheet_fname

		# Select lanes from flow cell that are bar coded.
		try:
			# Lane specification - mechanism to create sample sheets
			# for a subset of lanes on the flow cell.
			lane_specification=','.join(map(str,lanes))
			query = """select flowcell.barcode,
	flowcellchannel.number,
	sample.number,
	genomebuild.genomebuildname,
	sample.barcodesequence,
	appuser.firstname, 
	appuser.lastname,
	request.number,
	sample.barcodesequenceb,
	pipelineprotocol.protocol
from flowcell
	join flowcellchannel on flowcellchannel.idflowcell=flowcell.idflowcell
	join sequencelane on sequencelane.idflowcellchannel = 
		flowcellchannel.idflowcellchannel
	join sample on sequencelane.idsample = sample.idsample
	left outer join genomebuild on sequencelane.idgenomebuildalignto = 
		genomebuild.idgenomebuild
	join request on sequencelane.idrequest = request.idrequest
	join appuser on request.idappuser = appuser.idappuser
	left outer join pipelineprotocol on flowcellchannel.idpipelineprotocol = pipelineprotocol.idpipelineprotocol
where flowcellchannel.filename = '%s'
and flowcellchannel.number in (%s)
order by flowcellchannel.number, sample.number;""" % (self.id,lane_specification)
			c=self.Query(query)
		except pymssql.OperationError:
			print(query)
			raise
		results = c.fetchall()
		# Open file
		ofs = open(samplesheet_fname,'w')
		WriteSamplesheetHeader(ofs)
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
			fname=rec[5]
			lname=rec[6]
			rnum=rec[7]
			sample_barcode_b=rec[8]
			pipelineprotocol=rec[9]
			if sample_barcode_b is not None:
				sample_barcode = sample_barcode.strip()+'-'+sample_barcode_b.strip()
			elif type(sample_barcode) == type(''):
				sample_barcode = sample_barcode.strip()
			if pipelineprotocol is None:
				pipelineprotocol="None"
			# If lane has a single sample, don't write its
			# barcode. This turns off bar code processing for the 
			# lane for reasons described above.
			if sample_count[lane] == 1:
				sample_barcode = ''
			self.WriteSampleSheetRow( flowcell_barcode, lane, sampleid, genome, sample_barcode, fname, lname, rnum, pipelineprotocol, ofs )

		# Check if this is a rapid run. This will be the case
		# if the flow cell is registered as a single lane, but there
		# will be two lanes of data for it.
		#if len(sample_count.keys())==1 and \
		#	int(sample_count.keys()[0]) == 1 and \
		#	os.path.exists(os.path.join(self.dirname,'Data','Intensities','L002')):
		#	# Its a rapid run.
		#	self.log(["Run", self.id, "is a rapid run. Duplicating rows for lane 2 in sample sheet",samplesheet_fname,"."])
		#	for rec in results:
		#		# Write row.
		#		flowcell_barcode = rec[0]
		#		lane = rec[1]
		#		sampleid=rec[2]
		#		genome=rec[3]
		#		sample_barcode=rec[4]
		#		sample_barcode_b=rec[8]
		#		if sample_barcode_b is not None:
		#			sample_barcode = sample_barcode.strip()+'-'+sample_barcode_b.strip()
		#		elif type(sample_barcode) == type(''):
		#			sample_barcode = sample_barcode.strip()
		#		fname=rec[5]
		#		lname=rec[6]
		#		rnum=rec[7]
		#		# If lane has a single sample, don't write its
		#		# barcode.
		#		if sample_count[lane] == 1:
		#			sample_barcode = ''
		#		self.WriteSampleSheetRow( flowcell_barcode, 2, sampleid, genome, sample_barcode, fname, lname, rnum, ofs )

		# Close file.
		ofs.close()
		# Return file name.
		return samplesheet_fname

	def WriteSampleSheetRow( self, flowcell_barcode, lane, sampleid, genome, sample_barcode, fname, lname, rnum, pipelineprotocol, ofs ):
		"""
		Writes one row of a sample sheet file.
		Added pipeline protocol. Write this into the Recipe column.
		2017-07-10. BAM.
		"""
		row = [ flowcell_barcode,str(lane),sampleid,EraseCommas(genome) or 'None',sample_barcode or '',EraseCommas(fname+' '+lname),'N',pipelineprotocol,'Sandy',rnum ]
		ofs.write(','.join(row)+'\n')

	def InitializeCoreFacility( self ):
		"""
		InitializeCoreFacility identifies which core facility
		started this run.
		"""
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

	def IsProcessed(self):
		"""
		Returns True if sequencing run has already been processed. This
		is determined by finding the run in the sqlite database.
		"""
		return self.id in (State.complete,State.error)
# State machine. This defines which method to call depending on the
# Run's current state, and the Run's next state depending on whether
# the method returns True or False.
# The format of these entries is:
# current_state : ( method_to_call, state_if_True, state_if_False )
# By convention, the method names referenced in the state machine begin
# with Sm. This just makes it easier to find them in the code.

statemachine={
	State.new.value:( Run.SmCheckRegisteredVerbose,
		State.wait_for_data,
		State.wait_for_registration ),
	State.wait_for_registration.value:( Run.SmCheckRegisteredSilent,
		State.wait_for_data,
		State.wait_for_registration ),
	State.wait_for_data.value:( Run.SmCheckTransferComplete,
		State.pipeline,
		State.wait_for_data ),
	State.pipeline.value:( Run.SmRunPipeline,
		State.complete,
		State.error ),
	State.reprocess.value:( Run.SmReprocess,
		State.pipeline,
		State.error ),
	#State.complete.value:(method,t_state,f_state),
	#State.error.value:(method,t_state,f_state),
}

