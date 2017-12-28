import os
import xml.dom.minidom
import pipelineparams as params

from run import Run

class HiSeqRun(Run):
	
	def __init__( self, id, full_path_of_run_dir, state=None ):
		Run.__init__(self,id, full_path_of_run_dir, state)
		self.DetermineRunType()
	
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
			raise Exception("Expected 1 or 2 data reads in RunInfo.xml. Found %d data reads. Directory %s." % (self.num_data_reads,self.dirname ))

