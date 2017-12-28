#!/usr/bin/env python3
# run_illumina_pipeline.py - script to run low-level processing of illumina
# sequencer data.

import sys
from pipelinemgr import application

def main():
	try:
		app=application.Application()
		app.Go()
	except:
		sys.stderr.write("Exception!\n")
		raise


if __name__ == "__main__":
	main()
