"""
indexreads.py - definitions for IndexReads classes that define different 
numbers of index reads.
"""
class IndexReads:
	def NumIndexReads(self):
		return self.index_reads
	
	def IndexMaskComponents(self):
		return self.index_mask

class DualIndex(IndexReads):
	def __init__(self):
		self.index_reads=2
		self.run_type="Dual Index"
		if self.barcode_a_len > 0:
			barcode_a_mask = "I%dn*" % self.barcode_a_len
		else:
			barcode_a_mask = "n*"
		if self.barcode_b_len > 0:
			barcode_b_mask = "I%dn*" % self.barcode_b_len
		else:
			barcode_b_mask = "n*"
		self.index_mask=(barcode_a_mask,barcode_b_mask)

class SingleIndex(IndexReads):
	def __init__(self):
		self.index_reads=1
		self.run_type="Single Index"
		if self.barcode_a_len > 0:
			barcode_a_mask = "I%dn*" % self.barcode_a_len
		else:
			barcode_a_mask = "n*"
		self.index_mask=(barcode_a_mask,"")
