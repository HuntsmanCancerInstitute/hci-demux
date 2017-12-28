
import re
import os

from .logger import Logger
from . import pipelineparams

def ListRunDirectories():
	# Use this pattern to find both Hiseq and Miseq runs:
	pattern = re.compile("[0-9]*_[A-Z0-9]*_[0-9]*_[A-Z0-9-]*$")

	known_run_ids=[]
	new_runs=[]

	for subdir in os.listdir(pipelineparams.runs_dir_root):
		run_full_path=os.path.join(pipelineparams.runs_dir_root,subdir)
		if not os.path.isdir(run_full_path):
			# Skip non-directory files.
			continue
		if pattern.match(subdir):
			#self.Log(["Checking subdirectory",subdir])
			# Found a runs folder.
			if subdir not in known_run_ids:
				# Discovered a new run folder. 
				new_runs.append((subdir,run_full_path))
			else:
				Logger().Log(["Subdirectory",subdir,"already known."])
	new_runs.sort()
	return new_runs

def EraseCommas(s):
	"""Removes all commas in string s."""
	if s:
		return ''.join(s.split(','))
	else:
		return s

def WriteSamplesheetHeader(ofs):
	# Write header.
	header = ['FCID','Lane','SampleID','SampleRef','Index','Description','Control','Recipe','Operator','SampleProject']
	ofs.write(','.join(header)+'\n')
