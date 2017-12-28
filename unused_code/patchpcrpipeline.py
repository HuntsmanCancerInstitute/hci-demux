import xml.dom.minidom
import os
import re

from states import States
from runmgr import RunMgr
from simultaneousjobrunner import SimultaneousJobRunner

from patchpcr import PatchPcr

class PatchPcrPipeline(RunMgr,PatchPcr):
	def __init__(self, run_directories):
		self.type_of_runs = "PatchPcr"
		self.transition = {
			# curr_state, function, success_state, fail_state
			States.new: (PatchPcr.CheckRegisteredVerbose,
				States.wait_for_data,
				States.check_registered ),
			States.check_registered: (PatchPcr.CheckRegisteredSilent,
				States.wait_for_data,
				States.check_registered),
			States.wait_for_data: (PatchPcr.CheckTransferComplete,
				States.check_if_miseq,
				States.wait_for_data),
			States.check_if_miseq: (PatchPcr.CheckIfMiseq,
				States.check_if_patchpcr,
				States.error_detected),
			States.check_if_patchpcr: (PatchPcr.CheckIfPatchPcr,
				States.patchpcr_sample_sheet,
				States.error_detected),
			States.patchpcr_sample_sheet: (PatchPcr.PatchPcrSampleSheet,
				States.patchpcr_demultiplex,
				States.error_detected),
			States.patchpcr_demultiplex: (PatchPcr.PatchPcrDemultiplex,
				States.patchpcr_postprocess,
				States.error_detected),
			States.patchpcr_postprocess: (PatchPcr.PatchPcrPostprocess,
				States.patchpcr_distribute,
				States.error_detected),
			States.patchpcr_distribute: (PatchPcr.PatchPcrDistribute,
				States.Miseq_qc,
				States.error_detected),
			States.qc: (PatchPcr.Qc_18,
				States.archive,
				States.error_detected),
			States.archive: (PatchPcr.Archive,
				States.complete,
				States.error_detected),
			States.error_detected: (PatchPcr.NotifyError,
				States.error_notified,
				States.error_detected),
		}
		self.db_file = "pipelinemgr.db"

		#List of parent directories in which runs folders are created
		self.runs_to_process = run_directories

		#verify run directory exists
		for runs_dir in self.runs_to_process:
			if not os.path.exists(runs_dir):
				raise IOError, "Runs directory '%s does not exist." % runs_dir

		#Opens pipelinemgr database that tracks run processing. If it doesn't exist,
		# a new database is created
		if self.DbExists():
			self.DbOpen()
		else:
			self.DbCreate()

	def GetCorrectRuns(self):
		temp = self.GetAllRuns()
		for run in temp:
			run = PatchPcr(run[0],run[1],run[2])
			if run.CheckIfPatchPcr():
				self.runs_to_process.append(run)
