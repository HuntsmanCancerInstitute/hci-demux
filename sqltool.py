#!/usr/bin/env python

import hcidemux.gnomex
import hcidemux.pipelineparams as params
import sys
import re
import string
import unicodedata

def debugmsg(msg):
	sys.stderr.write("%s\n"%msg)

def CleanQuery( query ):
	query = re.sub("\n"," ",query)
	query = re.sub("\t"," ",query)
	#query = re.sub(" ","\n",query)
	return query

def ConvertToAscii(value):
	try:
		return(str(value))
	except UnicodeEncodeError:
		#sys.stderr.write("Encoding %s and ascii.\n" % value )
		return unicodedata.normalize('NFKD', value).encode('ascii','ignore')

def PrintRecords( recs ):
	#debugmsg(type(recs[0]))
	if type(recs[0]) == type({}):
		# Find column names.
		cols = []
		for item in recs[0].keys():
			if type(item)==type(''):
				cols.append(item)
		print string.join(cols,"\t")
		for rec in recs:
			vals=[]
			for col in cols:
				vals.append(str(rec[col]))
			print string.join(vals,"\t")
	else:
		for rec in recs:
			#print rec
			#print '\t'.join(map(str,list(rec)))
			print '\t'.join(map(ConvertToAscii,list(rec)))

def RunQuery( cursor, query, connection ):
	try:
		cursor.execute( query )
		# If its a query, return the results.
		if query[0:6] == 'select':
			recs = cursor.fetchall()
			if not recs:
				print 'No records selected.'
			else:
				PrintRecords( recs )
		elif query[0:6] == 'update':
			connection.commit()
			print '%d rows updated.' % cursor.rowcount
		else:
			print "Unknown type of database statement."
	except:
		print query
		raise

def ReadQuery():
	"""Reads a query from stdin. Reads until an empty line is found."""
	query = ""
	rec = sys.stdin.readline()
	while rec:
		#rec = rec.strip()
		if not rec.strip():
			break
		if query:
			query += ' '
		query += rec
		rec = sys.stdin.readline()
	return query

def RunQueries( cursor, connection ):
	while 1:
		query = ReadQuery()
		query = CleanQuery( query )
		if query:
			RunQuery( cursor, query, connection )
		else:
			break

def main():

	# Process any command-line arguments.
	debug=False
	asdict=False
	for arg in sys.argv[1:]:
		if arg == "--debug":
			debug=True
		elif arg == "--dict":
			asdict=True
	# Connect to db.
	g = hcidemux.gnomex.GNomExConnection()
	#connection = g.GnConnect()
	connection = g.GnConnect(params.db_user,params.db_password,asdict=False)
	cursor = connection.cursor()
	if debug:
		cursor._source.debug_queries = 1
		
	try:
		RunQueries( cursor, connection )
	except KeyboardInterrupt:
		sys.exit( 0 )

if __name__ == "__main__":
	main()
