

class Lane:
	"""Lane class represents one lane on a flow cell in a sequencing
	run. Lane has one method per state in state machine. Lane methods
	will be called by methods of Run. Special bioinformatics handling of 
	lanes can be implemented in classes derived from Lane.
	
	Examples of special handling are:
		save .cif files
		preserve all the bad barcode sequences
	"""

	def __init__( self, lane_number ):
		self.lane_number = lane_number

	def CheckSingleIndexLength( self ):
		print "Lane %d CheckSingleIndexLength." % self.lane_number

	def SimpleBclConvert( self ):
		print "Lane %d SimpleBclConvert." % self.lane_number

	def ComplexBclConvert( self ):
		print "Lane %d ComplexBclConvert." % self.lane_number

	def MergeBclConversions( self ):
		print "Lane %d MergeBclConversions." % self.lane_number

	def DistributeDemultiplexed( self ):
		print "Lane %d DistributeDemultiplexed." % self.lane_number

	def Cleanup( self ):
		print "Lane %d Cleanup." % self.lane_number

	def Qc( self ):
		print "Lane %d Qc." % self.lane_number

	def Archive( self ):
		print "Lane %d Archive." % self.lane_number

	def Complete( self ):
		print "Lane %d Complete." % self.lane_number

	def ErrorDetected( self ):
		print "Lane %d ErrorDetected." % self.lane_number

	def ErrorNotified( self ):
		print "Lane %d ErrorNotified." % self.lane_number

	def Reprocess( self ):
		print "Lane %d Reprocess." % self.lane_number

	state_method = {
		'e_check_single_index_length': CheckSingleIndexLength,
		'f_simple_bcl_convert': SimpleBclConvert,
		'g_complex_bcl_convert': ComplexBclConvert,
		'h_merge_bcl_conversions': MergeBclConversions,
		'i_distribute_demultiplexed': DistributeDemultiplexed,
		'j_cleanup': Cleanup,
		'k_qc': Qc,
		'l_archive': Archive,
		'm_complete': Complete,
		'n_error_detected': ErrorDetected,
		'o_error_notified': ErrorNotified,
		'p_reprocess': Reprocess,
	}

def test():
	lanes = []
	for lanenum in range(1,9):
		lanes.append(Lane(lanenum))
	
	states = Lane.state_method.keys()
	states.sort()
	for state in states:
		for lane in lanes:
			Lane.state_method[state]( lane )

if __name__ == "__main__":
	test()
