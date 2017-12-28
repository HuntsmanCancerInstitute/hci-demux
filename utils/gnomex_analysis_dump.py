#!/usr/bin/env python
# script to extract and dump analysis information from the GNomEx database for  
# a specific investigator, useful to combine with exported files from the Repository 
# before we do a big cleanup of their stuff
# writes out a tab delimited text file for importing into a spreadsheet

import sys
from hcidemux import gnomex
from hcidemux import pipelineparams as params

TAB="\t"


class ReportProducer():
	'''Produces a daily report of recently-sequenced samples.'''

	db_user=params.db_user
	db_password=params.db_password
	# Updated query to use codeapplication from the 
	# seqlibprotocolapplication table, since this is the one recorded
	# by the sequencing lab rather than the one requested by the user
	# (which may be incorrect).
	# Need to exclude the MiSeq runs, and select only samples from HiSeq 
	# runs.
	query="""select lab.firstname labfirst, 
		lab.lastname lablast, 
		appuser.firstname userfirst, 
		appuser.lastname userlast, 
		AnalysisGroup.name groupname, 
		Analysis.number analnumber, 
		Analysis.name analname, 
		organism.organism, 
		genomebuild.genomebuildname 
		from lab 
		join Analysis on Analysis.idLab = lab.idlab 
		join appuser on appuser.idappuser = Analysis.idappuser 
		join AnalysisGroupItem on AnalysisGroupItem.idAnalysis = Analysis.idAnalysis 
		left join AnalysisGroup on AnalysisGroupItem.idAnalysisGroup = AnalysisGroup.idAnalysisGroup 
		join organism on organism.idorganism = Analysis.idorganism 
		join AnalysisGenomeBuild on AnalysisGenomeBuild.idAnalysis = Analysis.idAnalysis 
		left join genomebuild on AnalysisGenomeBuild.idgenomebuild = genomebuild.idgenomebuild 
		where lab.firstname = '%s' and lab.lastname = '%s' 
		order by analnumber;"""

	#columns = ['Lab', 'Requester', 'Request', 'Project', 'Sample', 'Organism/Genome', 'Lane']
	def __init__( self ):
		# Connect to GNomEx
		g = gnomex.GNomExConnection()
		self.connection = g.GnConnect(ReportProducer.db_user,ReportProducer.db_password)
		self.cursor = self.connection.cursor()
		
	def RunReport(self,RequestFirstName,RequestLastName):
		self.cursor.execute(ReportProducer.query % (RequestFirstName, RequestLastName) )
		results = self.cursor.fetchall()
		# Document header.
		# Table header.
		columns = ['Lab', 'Owner', 'Group', 'Analysis', 'Analysis Name', 'Organism','Genome']
		doc = TAB.join( columns ) + '\n'

		# Table contents.
		for row in results:
			(labfirst, lablast, userfirst, userlast, groupname, analnumber, analname, organism, genomebuildname) = row
			if labfirst:
				labname = ' '.join( [labfirst,lablast,'Lab'])
			else:
				labname = lablast
			username = ' '.join([userfirst,userlast])
			if not genomebuildname:
				genomebuildname = 'None'
			doc+=TAB.join([labname,username,groupname,analnumber,analname,organism,genomebuildname]) + '\n'

		# Send the report.
		self.WriteReport( doc, RequestFirstName, RequestLastName )

	def WriteReport( self, report, RequestFirstName, RequestLastName ):
		# Create a file name.
		outfile = '%s_%s_Analysis.txt' % ( RequestFirstName, RequestLastName )
		ofs = open( outfile, 'w' )
		ofs.write( report )
		ofs.close()
		sys.stderr.write(" Report written to %s.\n" % ( outfile ))

def Usage():
	sys.stderr.write("Use: %s firstName lastName \nwhere the name is the name of the investigator.\n" % sys.argv[0])

def main():
	if len(sys.argv) == 3:
		RequestFirstName = sys.argv[1]
		RequestLastName = sys.argv[2]
		ReportProducer().RunReport(RequestFirstName,RequestLastName)
	else:
		Usage()
		sys.exit(0)

if __name__ == "__main__":
	main()
