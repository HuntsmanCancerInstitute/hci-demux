import os
import re
import sys

from states import States
from runmgr import RunMgr
from simultaneousjobrunner import SimultaneousJobRunner

from kappapcr import KappaPcr

class KappaPcrPipeline(RunMgr,KappaPcr):
	def __init__(self, run_directories):
		self.type_of_runs = "KappaPcr"
		self.transition = {
			# curr_state, function, success_state, fail_state
			States.new: (KappaPcr.CheckRegisteredVerbose,
				States.wait_for_data,
				States.check_registered ),
			States.check_registered: (KappaPcr.CheckRegisteredSilent,
				States.wait_for_data,
				States.check_registered),
			States.wait_for_data: (KappaPcr.CheckTransferComplete,
				States.check_if_miseq,
				States.wait_for_data),
			States.check_if_miseq: (KappaPcr.CheckIfMiseq,
				States.check_if_kappapcr,
				States.error_detected),
			States.check_if_kappapcr: (KappaPcr.CheckIfKappaPcr,
				States.kappapcr_sample_sheet,
				States.error_detected),
			States.kappapcr_sample_sheet: (KappaPcr.KappaPcrSampleSheet,
				States.kappapcr_demultiplex,
				States.error_detected),
			States.kappapcr_demultiplex: (KappaPcr.KappaPcrDemultiplex,
				States.kappapcr_postprocess,
				States.error_detected),
			States.kappapcr_postprocess: (KappaPcr.KappaPcrPostprocess,
				States.kappapcr_distribute,
				States.error_detected),
			States.kappapcr_distribute: (KappaPcr.KappaPcrDistribute,
				States.qc,
				States.error_detected),
			States.qc: (KappaPcr.Qc_18,
				States.archive,
				States.error_detected),
			States.archive: (KappaPcr.Archive,
				States.complete,
				States.error_detected),
			States.error_detected: (KappaPcr.NotifyError,
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
			run = KappaPcr(run[0],run[1],run[2])
			if run.CheckIfKappaPcr():
				self.runs_to_process.append(run)
