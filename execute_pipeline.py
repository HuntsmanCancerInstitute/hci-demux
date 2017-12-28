#!/usr/bin/env python
# execute_pipeline.py - top level code to run hci_demux pipelines.
import sys
import os
import xml.dom.minidom

try:
	import sysv_ipc
except ImportError:
	Logger().Log( "%s requires module sysv_ipc. See http://semanchuk.com/philip/sysv_ipc/." % sys.argv[0] )
	sys.exit(1)

from hcidemux.runmgr import RunMgr
from hcidemux.logger import Logger
import hcidemux.pipelineparams as params
import hcidemux.gnomex as gnomex

# Run classes:
from hcidemux.patchpcrrun import PatchPcrRun
from hcidemux.kappapcrrun import KappaPcrRun
from hcidemux.singleendhiseqrun import SEHiSeqRun
from hcidemux.pairedendhiseqrun import PEHiSeqRun

def CountRunDataReads(run_full_path):
	"""
	Given the full path name of a run folder, this function returns the
	number of data reads for the sequencing run by parsing the RunInfo.xml
	file.
	"""
	# Parse the RunInfo.xml file.
	configfile = os.path.join(run_full_path, "RunInfo.xml")
	doc = xml.dom.minidom.parse(configfile)

	# Get the Read elements from the document.
	reads=doc.getElementsByTagName("Read")

	# Count the reads that are at least 25 bp long.
	num_data_reads = 0
	for read in reads:
		if int(read.getAttribute('NumCycles')) >= 25:
			num_data_reads += 1
	return num_data_reads

def IdentifyPipeline(run_id,run_full_path):
	"""
	IdentifyPipeline(run_id) - given id of a sequencing run (which is the
	directory name of the run folder, not full path) returns 
	the class that is to process the run, or None of the type of run
	isn't recognized.
	"""
	# Connect to GNomex.
	connection = gnomex.GNomExConnection().GnConnect(params.db_user,params.db_password)
	c = connection.cursor()
	# Select the sequencing application name(s) in use on the flow cell.
	query="""select distinct application.application
		from application 
		join seqlibprotocolapplication 
		on application.codeapplication = seqlibprotocolapplication.codeapplication
		join seqlibprotocol on seqlibprotocol.idseqlibprotocol = seqlibprotocolapplication.idseqlibprotocol
		join sample on sample.idseqlibprotocol = seqlibprotocol.idseqlibprotocol
		join sequencelane on sequencelane.idsample = sample.idsample
		join flowcellchannel on sequencelane.idflowcellchannel = flowcellchannel.idflowcellchannel
		where flowcellchannel.filename = '%s';"""% run_id
	c.execute(query)
	recs=c.fetchall()
	# recs is a list of tuples each containing a unicode string. Convert
	# this into a single pipe delimited string.
	applications = '|'.join(map( lambda x: str(x[0]), recs ))

	Logger().Log("Run %s seq applications: %s." % ( run_id, applications ))
	num_reads = CountRunDataReads(run_full_path)
	Logger().Log("Run %s has %d data reads." % (run_id,num_reads))
	if applications.find("Patch PCR") > -1:
		return PatchPcrRun
	elif applications.find("Kappa PCR") > -1:
		return KappaPcrRun
	elif num_reads == 2:
		# This is a standard paired-end HiSeq run.
		return PEHiSeqRun
	elif num_reads == 1:
		# This is a standard single-end HiSeq run.
		return SEHiSeqRun
	return None

def execute(pipeline_class):

	try:
		mgr = pipeline_class( params.root_directories )
		print(params.root_directories)
		mgr.Discover()
		mgr.GetCorrectRuns()
		Logger().Log("Grabbing all %s runs" % mgr.type_of_runs)
		mgr.ProcessRuns()
	except IOError, message:
		print("IOERROR")
		print(params.root_directories)
		mgr.Discover()
		mgr.GetCorrectRuns()
		Logger().Log("%s - a runs directory in %s doesn't exist." % ( sys.argv[0], params.root_directories ))
	except:
		# Clean up after unforseen errors.
		Logger().Log("%s - unexpected exception caught." % sys.argv[0])
		sem.release()
		sem.remove()
		raise


def main():
	"""
	main() function for hci_demux package.
	"""
	# Set up semaphore to guarantee exclusivity.
	keyval = 43
	flags = sysv_ipc.IPC_CREX
	initial = 1
	try:
		sem = sysv_ipc.Semaphore(key=keyval, flags=flags, initial_value=initial)
		sem.acquire(timeout=0)
	except sysv_ipc.ExistentialError:
		Logger().Log( "%s - program already running." % sys.argv[0] )
		sys.exit(0)

	try:
		# Create a generic RunMgr object.
		run_mgr = RunMgr()

		# Discover new sequencing runs and add them to the sqlite 
		# database.
		run_mgr.Discover(params.root_directories)

		# Process active runs, ie runs not finished or in error state.
		for (run_id,run_full_path,run_state) in run_mgr.GetActiveRuns():
			# Identify the class to handle the run and the class
			# of the run itself.
			run_class = IdentifyPipeline(run_id,run_full_path)
			# Create an object of that specific pipeline class.
			if run_class is None:
				# We could not find a pipeline for the run.
				Logger().Log("No pipeline class found for run %s." % run_full_path)
			else:
				Logger().Log("Will process run %s with run class %s." % ( run_id, run_class ) )
				run_object = run_class(run_id,run_full_path,run_state)
				# Process the run.
				run_mgr.AddRun( run_object )
		run_mgr.ProcessRuns()
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
