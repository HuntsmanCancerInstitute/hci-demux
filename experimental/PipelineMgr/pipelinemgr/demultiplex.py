"""
demultiplex.py - definitions for classes that define different demultiplexing
schemes.
"""
class Demultiplex:
	def BclConfigMismatches(self):
		return ["--barcode-mismatches",str(self.num_mismatches)]

	def NumMismatches(self):
		return self.num_mismatches

	def DemultiplexMethod(self):
		return self.demultiplex_method

class SingleMismatch(Demultiplex):
	def __init__(self):
		self.num_mismatches=1
		self.demultiplex_method = "Single Mismatch"
	
class NoMismatch(Demultiplex):
	def __init__(self):
		self.num_mismatches=0
		self.demultiplex_method = "No Mismatch"
	
class MolecularBarcode(Demultiplex):
	def __init__(self):
		self.num_mismatches=1
		self.demultiplex_method = "Molecular Barcode"
		# This class redefines the index mask because one of the
		# index reads becomes a data read.
		if self.barcode_a_len > 0:
			barcode_a_mask = "I%dn*" % self.barcode_a_len
		else:
			barcode_a_mask = "n*"
		barcode_b_mask = "Y*"
		self.index_mask=(barcode_a_mask,barcode_b_mask)

class NoDemultiplex(Demultiplex):
	def __init__(self):
		self.num_mismatches=0
		self.demultiplex_method = "No Demultiplex"

class UnknownDemultiplexMethodException(Exception):
	pass

# Map each demultiplex class to a name. These names are in the 
# PiplineProtocols table in GNomEx.
demux_methods= {
	'Single base mismatch': SingleMismatch,
	'Zero base mismatch': NoMismatch,
	'No demultiplexing': NoDemultiplex,
	'Second barcode random N-mer': MolecularBarcode,
}
