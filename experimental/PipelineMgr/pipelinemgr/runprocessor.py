#!/usr/bin/env python3
# RunProcessor - code to automate processing Illumina sequencing runs.

import sys
import socket
import time
import re
import os
import subprocess

from . import gnomex
from . import pipelineparams
from .logger import Logger
from .emailer import Emailer

from .datareads import PairedEnd, SingleEnd
from .indexreads import DualIndex, SingleIndex
from .instrument import HiSeq, MiSeq

from .instrument import UnknownInstrumentException
from .demultiplex import UnknownDemultiplexMethodException
from .demultiplex import demux_methods
from .util import EraseCommas

# Class that represents one sample in the sequencing run:

class SeqSample(object):
	"""
	SeqSample represents one sample in a sequencing run. Together, all the
	SeqSample objects for a run contain all the information needed to create
	the sample sheet for the run.
	"""
	def __init__(self, flowcell_barcode, lane, sample_id, genome, sample_barcode, first_name, last_name, request_number ):
		self.flowcell_barcode=flowcell_barcode
		self.lane=lane
		self.sample_id=sample_id
		self.genome=genome
		self.sample_barcode=sample_barcode
		self.first_name=first_name
		self.last_name=last_name
		self.request_number=request_number

	def WriteSampleSheetRow( self, ofs ):
		# Regular row. Write it to the file.
		row = [ self.flowcell_barcode,repr(self.lane),self.sample_id,EraseCommas(self.genome) or 'None',self.sample_barcode or '',EraseCommas(self.first_name+' '+self.last_name),'N','RI','Sandy',self.request_number ]
		ofs.write(','.join(row)+'\n')

# Class providing the run processing capability:

class RunProcessor(Logger,Emailer):
	"""
	A RunProcessor object handles BCL conversion, demultiplexing, 
	and data distribution for one or more lanes of a flow cell that have
	identical characteristics, ie # of data reads, # of bar code reads,
	bar code lengths, demultiplexing strategy, etc.
	"""

	def __init__(self,i,run_name,run_dir,num_data_reads,num_index_reads,sample_sheet_name,barcode_a_len,barcode_b_len,samples=[]):
		self.i=i
		self.id=run_name
		self.run_directory=run_dir
		#self.gnomex_connection = GNomEx.GNomExConnection().GnConnect(pipelineparams.db_user,pipelineparams.db_password)
		self.samples=samples[:]
		self.lanes=[]
		self.num_data_reads=num_data_reads
		self.num_index_reads=num_index_reads
		self.sample_sheet_name=sample_sheet_name
		self.barcode_a_len=barcode_a_len
		self.barcode_b_len=barcode_b_len

		#self.InitializeCoreFacility()
	
	def AddLane(self,lane):
		self.lanes.append(lane)

	def AddSample(self,lane,sampleid,index,project):
		pass

	def FromAddr( self ):
		"""
		Returns from address for emailing purposes.
		"""
		# Unix account under which PipelineMgr.py runs.
		username = os.environ['LOGNAME']
		return username + "@" + socket.getfqdn()

	def Describe(self):
		print("Run directory:", self.run_directory)
		self.PrintInstrument()
		self.PrintRunType()
		self.PrintDemultiplexingScheme()
		print("Sequencer finished:", self.SequencerFinished())
		print()

	#def Query( self, query ):
	#	"""
	#	Runs a SQL query against GNomEx, and retuns a cursor from
	#	which the results can be retrieved.
	#	"""
	#	c = self.gnomex_connection.cursor()
	#	c.execute(query)
	#	return c

	def ProcessRun(self):
		"""
		Main method to call to process a run. This method calls all
		the steps necessary, many of which will be provided by the other
		classes.
		Argument i is the number of RunProcessor, assigned when
		"""
		# Make sample sheet. - RunProcessor class.
		# Convert BCL files. - Instrument class.
		if not self.ConvertBclFiles():
			self.Log("Run %s: BCL file conversion problem." % self.run_directory)
			return False

		# Merge BCL file batches.
		# Rename data files.
		# Distribute data files.
		# Generate QC report.
		# Notify that run is complete.

		self.Log("Run %s processing completed." % self.run_directory)
		return True
	
	def BclConfigInputDir(self):
		return ["--input-dir",os.path.join(self.run_directory,"Data","Intensities","BaseCalls")]
	
	def BclConfigOutputDir(self):
		self.lanes.sort()
		#outdir="_".join(["Unaligned"]+map(str,self.lanes))
		outdir="Unaligned_%d" % self.i
		outdir_full_path=os.path.join(self.run_directory,outdir)
		return ["--output-dir",outdir_full_path]

	def BclConfigUseBasesMask(self):
		(d1,d2) = self.DataMaskComponents()
		(i1,i2) = self.IndexMaskComponents()
		bases_mask=','.join(filter(lambda x:len(x)>0, [d1,i1,i2,d2]))
		return ["--use-bases-mask",bases_mask]

	def BclConfigCommonFlags(self):
		return ["--ignore-missing-bcls"]

	def BclConfigSampleSheet(self):
		return ["--sample-sheet",os.path.join(self.run_directory,"created_samplesheet_%d.csv" % self.i)]

	def ConvertBclFiles(self):
		"""
		ConvertBclFiles call bcl2fastq to perform the BCL to Fastq
		conversion.
		"""
		application_name="bcl2fastq"
		child_dir=self.run_directory
		child_stdout=open(os.path.join(child_dir,"bcltofastq_%d.out"%self.i),'w')
		child_stderr=open(os.path.join(child_dir,"bcltofastq_%d.err"%self.i),'w')
		child_args=[os.path.join( pipelineparams.bcl2fastq2_dir ,application_name) ] \
			+ self.BclConfigInputDir() \
			+ self.BclConfigOutputDir() \
			+ self.BclConfigSampleSheet() \
			+ self.BclConfigMismatches() \
			+ self.BclConfigUseBasesMask() \
			+ self.BclConfigCommonFlags()

		self.Log("Running configBclToFastq.")
		self.Log("\n\t".join(child_args))
		p = subprocess.Popen(args=child_args,cwd=child_dir,stdout=child_stdout,stderr=child_stderr)
		retval = p.wait()
		child_stdout.close()
		child_stderr.close()
		self.Log("%s returned %d." % (application_name,retval))
		return retval == 0


#def GetRunSamples(run_id):
#	"""
#	Queries GNomEx database for information about all samples in
#	a sequencing run. Returns a list of SeqSample objects.
#	"""
#	seqsamples=[]
#	gnomex_connection = GNomEx.GNomExConnection().GnConnect(pipelineparams.db_user,pipelineparams.db_password)
#
#	# Select lanes from flow cell that are bar coded.
#	c = gnomex_connection.cursor()
#	try:
#		query = """select flowcell.barcode,
#flowcellchannel.number,
#sample.number,
#genomebuild.genomebuildname,
#sample.barcodesequence,
#appuser.firstname, 
#appuser.lastname,
#request.number,
#sample.barcodesequenceb
#from flowcell
#join flowcellchannel on flowcellchannel.idflowcell=flowcell.idflowcell
#join sequencelane on sequencelane.idflowcellchannel = 
#	flowcellchannel.idflowcellchannel
#join sample on sequencelane.idsample = sample.idsample
#left outer join genomebuild on sequencelane.idgenomebuildalignto = 
#	genomebuild.idgenomebuild
#join request on sequencelane.idrequest = request.idrequest
#join appuser on request.idappuser = appuser.idappuser
#where flowcellchannel.filename = '%s'
#order by flowcellchannel.number, sample.number;""" % run_id
#		c.execute(query)
#	except pymssql.OperationError:
#		print(query)
#		raise
#	results = c.fetchall()
#	for rec in results:
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
#		seqsamples.append(SeqSample(flowcell_barcode,lane,sampleid,genome,sample_barcode,fname,lname,rnum))
#	c.close()
#	return seqsamples

def CreateRunProcessor(i,run_full_path,instrument,num_data_reads,num_index_reads,barcode_a_len,barcode_b_len,demultiplex_method,sample_sheet_name):
	"""
	Given the instrument type, number of data reads, barcode a length,
	barcode b length, and demultiplexing method, returns a RunProcessor
	object to handle one or more lanes.
	After the object is created, the code must add lanes to the object
	before it can process anything.
	"""
	Logger().Log("Creating a run processor for %s." % run_full_path)

	assert(num_data_reads in [1,2])
	assert(instrument.lower() in ['hiseq','miseq'])

	# Get a list of samples on the run.
	# Do this later - once all the lanes have been added to the 
	# run processor.
	#samples=GetRunSamples(run_name)

	# Determine instrument for sequencing run. Using the LaneCount attribute
	# of the FlowcellLayout tag in the RunInfo.xml file for this.
	if instrument.lower() == 'miseq':
		InstrumentClass=MiSeq
	elif instrument.lower() == 'hiseq':
		InstrumentClass=HiSeq

	# Determine the number of Datareads class.
	if num_data_reads == 1:
		DatareadsClass=SingleEnd
	elif num_data_reads == 2:
		DatareadsClass=PairedEnd
	
	# Determine the number of Indexreads class.
	if num_index_reads == 1:
		IndexreadsClass = SingleIndex
	elif num_index_reads == 2:
		IndexreadsClass = DualIndex

	# Determine the barcoding scheme.
	try:
		DemuxClass=demux_methods[demultiplex_method]
	except KeyError:
		raise UnknownDemultiplexMethodException("No demultiplexing class for method '%s' found." % demultiplex_method)

	# Make the RunProcessor class for the run. The DemuxClass must come 
	# after the IndexreadsClass because it may override some IndexreadsClass
	# values.
	run_name=os.path.split(run_full_path)[-1]
	objclass=type('runproc_'+run_name,(RunProcessor,InstrumentClass,DatareadsClass,IndexreadsClass,DemuxClass,),{})

	# Make an instance of the RunProcessor for this run.
	runproc_instance=objclass(i,run_name,run_full_path,num_data_reads,num_index_reads,sample_sheet_name,barcode_a_len,barcode_b_len)
	# For some reason we need to explicitly call the base class __init__ 
	# functions for all but the first base class.
	for baseclass in runproc_instance.__class__.__bases__[1:]:
		baseclass.__init__(runproc_instance)

	# Don't process the run yet. The code calling this function from
	# run.py needs to add some lanes to the runprocessor first.
	#runproc_instance.ProcessRun()
	return runproc_instance
