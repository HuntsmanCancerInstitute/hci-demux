#!/usr/bin/env python
# script to extract and dump sequencing request information from the GNomEx database for  
# a specific investigator, useful to combine with exported files from the Repository 
# before we do a big cleanup of their stuff
# writes out a tab delimited text file for importing into a spreadsheet

import sys
from hcidemux import gnomex
import hcidemux.pipelineparams as params

TAB="\t"


def HtmlTd( content ):
	return '<TD>%s</TD>' % content

def HtmlTh( content ):
	return '<TH>%s</TH>' % content

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
		request.number reqnum,
		project.name projname,
		sample.number samnum, 
		sample.name samname,
		organism.organism,
		genomebuild.genomebuildname,
		seqlibprotocolapplication.codeapplication,
		request.createdate,
		numbersequencingcycles.numbersequencingcycles,
		seqruntype.seqruntype
	from lab
	join request on request.idlab = lab.idlab
	join appuser on appuser.idappuser = request.idappuser
	join project on request.idproject = project.idproject
	join sample on sample.idrequest = request.idrequest
	join organism on sample.idorganism = organism.idorganism
	join sequencelane on sequencelane.idsample = sample.idsample
	left outer join genomebuild on sequencelane.idgenomebuildalignto = genomebuild.idgenomebuild
	join seqruntype on sequencelane.idseqruntype = seqruntype.idseqruntype
	join numbersequencingcycles on sequencelane.idnumbersequencingcycles = numbersequencingcycles.idnumbersequencingcycles
	join seqlibprotocol on sample.idseqlibprotocol = seqlibprotocol.idseqlibprotocol
	join seqlibprotocolapplication on seqlibprotocol.idseqlibprotocol = seqlibprotocolapplication.idseqlibprotocol
	join application on seqlibprotocolapplication.codeapplication = application.codeapplication
	where lab.firstname = '%s'
	and lab.lastname = '%s'
	order by userlast, samnum;"""

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
		columns = ['Lab', 'Requester', 'RequestNum', 'Project', 'SampleId', 'SampleName', 'Organism','Genome','Cycles','RunType','RequestYear']
		doc = TAB.join( columns ) + '\n'

		# Table contents.
		for row in results:
			(labfirst, lablast, userfirst, userlast, reqnum, projname, samnum, samname, organism, genomebuildname,applicationcode,reqdate,numcycles,seqruntype) = row
			if labfirst:
				labname = ' '.join( [labfirst,lablast,'Lab'])
			else:
				labname = lablast
			username = ' '.join([userfirst,userlast])
			sample = samnum + '<BR>'+samname
			if not genomebuildname:
				genomebuildname = 'None'
			doc+=TAB.join([labname,username,reqnum,projname,samnum,samname,organism,genomebuildname,`numcycles`,seqruntype,`reqdate.year`]) + '\n'

		# Send the report.
		self.WriteReport( doc, RequestFirstName, RequestLastName )

	def WriteReport( self, report, RequestFirstName, RequestLastName ):
		# Create a file name.
		outfile = '%s_%s_Request.txt' % ( RequestFirstName, RequestLastName )
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
