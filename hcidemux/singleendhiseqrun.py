from hiseqrun import HiSeqRun
from states import States

class SEHiSeqRun(HiSeqRun):
	def __init__( self, id, full_path_of_run_dir, state=None ):
		HiSeqRun.__init__(self,id, full_path_of_run_dir, state)
		self.type_of_runs = "Single End HiSeq"
		self.isPaired = False
		self.transition = {
			# curr_state, function, success_state, fail_state
			States.new: (HiSeqRun.CheckRegisteredVerbose,
				States.wait_for_data,
				States.check_registered ),
			States.check_registered: (HiSeqRun.CheckRegisteredSilent,
				States.wait_for_data,
				States.check_registered),
			States.wait_for_data: (HiSeqRun.CheckTransferComplete,
				States.make_sample_sheet,
				States.wait_for_data),
			States.make_sample_sheet: (HiSeqRun.MakeSampleSheet,
				States.check_single_index_length,
				States.error_detected),
			States.check_single_index_length: (HiSeqRun.CheckSingleIndexLength,
				States.simple_bcl_convert,
				States.complex_bcl_convert),
			States.simple_bcl_convert: (HiSeqRun.BclConvert,
				States.distribute_demultiplexed,
				States.error_detected),
			States.complex_bcl_convert: (HiSeqRun.BclConvertComplex,
				States.merge_bcl_conversions,
				States.error_detected),
			States.merge_bcl_conversions: (HiSeqRun.MergeBclConversions,
				States.distribute_demultiplexed,
				States.error_detected),
			States.distribute_demultiplexed: (HiSeqRun.DistributeDemultiplexed_18,
				States.cleanup,
				States.error_detected),
			States.cleanup: (HiSeqRun.Cleanup,
				States.qc,
				States.error_detected),
			States.qc: (HiSeqRun.Qc_18,
				States.archive,
				States.error_detected),
			States.archive: (HiSeqRun.Archive,
				States.complete,
				States.error_detected),
			States.error_detected: (HiSeqRun.NotifyError,
				States.error_notified,
				States.error_detected),
			States.reprocess: (HiSeqRun.Reprocess,
				States.make_sample_sheet,
				States.error_detected ),
		}

def main():
	se_run = SEHiSeqRun("bar","/foo/bar")

if __name__ == "__main__":
	main()
