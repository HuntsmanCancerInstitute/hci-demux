#!/usr/bin/env python

import pymssql
import re
import os

class GNomExConnection:

	def CleanQuery( self, query ):
		query = re.sub("\t"," ",query)
		return query

	def GnConnect(self,db_user,db_password,asdict=False):
		'''
		Returns a database connection. Next step is to call its 
		cursor() method, then cursor.execute(query) and cursor.fetchall()
		for example.
		'''
		
		try:
			host=os.environ['DBHOST']
		except KeyError:
			host='hci-db.hci.utah.edu:1433'
		database='GNomEx'
		return pymssql.connect(host=host,user=db_user,password=db_password,database=database,as_dict=asdict)
	
def test():
	g = GNomExConnection()
	connection = g.GnConnect(db_user="bmilash",db_password="123zxc",asdict=False)
	c = connection.cursor()
	c.execute("select count(*) n from AppUser")
	results = c.fetchall()
	try:
		print "There are %d registered GNomEx users." % results[0][0]
	except:
		print "results=",results
		raise

if __name__ == "__main__":
	test()
