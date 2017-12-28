"""
datareads.py - definitions for DataReads classes that define different 
numbers of data reads.
"""
class DataReads:
	def NumDataReads(self):
		return self.data_reads
	
	def DataMaskComponents(self):
		return self.data_mask

class PairedEnd(DataReads):
	def __init__(self):
		self.data_reads=2
		self.run_type="Paired End"
		self.data_mask=("Y*","Y*")

class SingleEnd(DataReads):
	def __init__(self):
		self.data_reads=1
		self.run_type="Single End"
		self.data_mask=("Y*","")
