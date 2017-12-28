"""
instrument.py - classes that define different Instrument classes.
"""

class Instrument:
	def InstrumentName(self):
		return self.instrument_name
	
	def BclConfigPositionsFormat(self):
		return ["--positions-format",self.positions_format]

class HiSeq(Instrument):
	def __init__(self):
		self.positions_format=".clocs"
		self.instrument_name="HiSeq"

class MiSeq(Instrument):
	def __init__(self):
		self.positions_format=".locs"
		self.instrument_name="MiSeq"

class UnknownInstrumentException(Exception):
	pass
