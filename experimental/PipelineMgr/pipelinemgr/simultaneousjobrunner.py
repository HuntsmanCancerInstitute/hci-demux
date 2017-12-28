#!/usr/bin/env python3
import os
import sys
import string
import subprocess

from .logger import Logger

class SimultaneousJobRunner(Logger):
	"""
	SimultaneousJobRunner runs a list of tasks using a fixed number
	of subprocesses. The number of simultaneous subprocesses is determined
	by the RunJobs method. New jobs are added to the list using the AddJob 
	method. This class is derived from Logger, which handles logging
	any messages. The main() function also works as a command-line
	utility to run a list of jobs read from a file or stdin.
	"""
	def __init__(self, verbose=False,save_output=False):
		#self.jobs = []
		#self.pids = {}
		#self.pipes = {}
		#self.output = {}
		self.verbose = verbose
		self.save_output = save_output
		#self.nextjob=0
		self.ClearJobs()
		#if self.verbose:
		#	self.log( "In SimultaneousJobRunner.__init__")

	def ClearJobs( self ):
		"""
		Initializes and clears out any previous batch of jobs.
		"""
		self.jobs = []
		self.pids = {}
		self.pipes = {}
		self.output = {}
		self.nextjob=0

	def AddJob( self, statement, output_key='' ):
		if self.verbose:
			self.log(["Adding job", statement])
		self.jobs.append( statement )

	def StartJob( self, verbose=False ):
		"""Starts the next job on the list."""
		if self.nextjob < len(self.jobs):
			job = self.jobs[ self.nextjob ]
			if self.verbose or verbose:
				self.log( "Starting job %d of %d, %s" % (self.nextjob,len(self.jobs),job))
			self.nextjob+=1

			# Start the subprocess.
			pipe=subprocess.Popen(job,shell=True)
			self.pipes[pipe.pid] = pipe
			# Save its process id.
			self.pids[pipe.pid] = job

	def RunJobs( self, numjobs, verbose=False ):
		"""Run up to numjobs at a time. As soon as one finishes, start another."""
		batchexitstatus=0

		# Start initial jobs.
		for i in range(0,numjobs):
			self.StartJob()

		# Wait for a job to finish, and start the next job.
		while self.pids:
			if self.verbose or verbose:
				pids = list(self.pids.keys())
				pids.sort()
				self.log( "Waiting on %d jobs: %s" % (len(self.pids), repr(pids)))
			(pid,jobexitstatus) = os.waitpid(-1,os.P_WAIT)
			jobname = self.pids.pop(pid)
			if self.verbose or verbose:
				self.log( "Pid %d (%s) finished, exitstatus=%d." % ( pid, jobname, jobexitstatus))
				if jobexitstatus != 0:
					batchexitstatus = jobexitstatus
			try:
				self.StartJob(verbose)
			except IndexError:
				pass

			# Confirm that the running pids are really running,
			# and start additional jobs as needed.
			for pid in list(self.pids.keys()):
				try:
					os.kill(pid,0)
				except OSError:
					# That pid is no longer running.
					self.pids.pop(pid)
					self.log( "Pid %d (%s) finished, exitstatus unknown." % ( pid, jobname ))
					try:
						self.output[pid] = self.pipes[pid].read()
					except AttributeError:
						pass
					self.StartJob()
		self.jobs = []
		self.nextjob=0
		self.log("All jobs complete. Job queue is empty.")
		return batchexitstatus == 0

def main():
	# Defaults: run 5 simultaneous jobs, reads jobs from stdin.
	numjobs=5
	ifs=sys.stdin
	verbose=False

	# Command-line argument processing.
	for i in range(1,len(sys.argv)):
		if sys.argv[i] == "-v":
			verbose=True
		if sys.argv[i] == "-n":
			try:
				numjobs = int(sys.argv[i+1])
			except IndexError:
				pass
		if os.path.exists(sys.argv[i]):
			ifs =open(sys.argv[i])
	
	jobs = list(map(string.strip,ifs.readlines()))
	if ifs != sys.stdin:
		ifs.close()
	
	Logger().log(["Verbose",verbose,"Numjobs",numjobs])
	s = SimultaneousJobRunner(verbose)
	for job in jobs:
		s.AddJob(job)
	if s.RunJobs( numjobs ):
		s.log("All jobs completed successfully.")
	else:
		s.log("One or more jobs finished with non-zero exit status.")

if __name__ == "__main__":
	main()
