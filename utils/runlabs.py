#!/usr/bin/env python
# runlabs.py - for a given run, list the labs represented on the flow cell.

import sys
import os
from hcidemux import gnomex
from hcidemux import pipelineparams as params

def ReportRunLabs( run_folder_name ):
	"""Given a run directory name, list the request numbers and the 
	requesting lab names found in the run's sample sheet file."""
	# Locate the run sample sheet.
	sample_sheet = os.path.join( run_folder_name, "Data","Intensities","BaseCalls","created_samplesheet.csv")
	#print "Sample sheet:", sample_sheet
	lanes = {}
	try:
		ifs=open(sample_sheet)
		# Collect list of request numbers for run.
		requests = set()
		# Read and discard header.
		ifs.readline()
		for rec in ifs:
			f = rec.strip().split(',')
			request = f[9]
			lane = f[1]
			requests.add(request)
			try:
				lanes[request].add(lane)
			except KeyError:
				lanes[request] = set()
				lanes[request].add(lane)
		ifs.close()
		#print requests
		# Select request number and lab name.
		request_string = ",".join(map(lambda x:"'"+x+"'",list(requests)))
		query = """select request.number,lab.firstname,lab.lastname
from request,lab
where request.idlab=lab.idlab
and request.number in (%s)"""  % request_string
		#print query
		g = gnomex.GNomExConnection()
		connection = g.GnConnect(params.db_user,params.db_password,asdict=True)
		c = connection.cursor()
		c.execute(query)
		results = c.fetchall()
		print
		print "Run folder:", run_folder_name
		for rec in results:
			reqnum = rec['number']
			firstname = rec['firstname']
			lastname = rec['lastname']
			if firstname:
				name = firstname + ' ' + lastname
			else:
				name = lastname
			lane_list = list(lanes[reqnum])
			lane_list.sort()
			lane = ','.join(lane_list)
			print '\t'.join([reqnum,name,lane])

	except IOError:
		sys.stderr.write("Sample sheet %s not found.\n" % sample_sheet)
		sys.exit(2)

def Usage():
	sys.stderr.write("Use: %s run_folder_name [ run_folder_name ... ]\n" % sys.argv[0] )

def main():
	# Check command-line arguments. If none, print the usage message.
	if len(sys.argv) < 2:
		Usage()
		sys.exit(1)
	
	for run_folder in sys.argv[1:]:
		ReportRunLabs( run_folder )

if __name__ == "__main__":
	main()
