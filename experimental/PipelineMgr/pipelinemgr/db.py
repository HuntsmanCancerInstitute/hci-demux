"""
db.py - database-related code. Utility functions for creating, opening, and
closing the database, as well as a base class for sequencing run objects
that need to save themselves to the database.
"""
import os
try:
	import sqlite3
except ImportError:
	import sqlite as sqlite3

from . import pipelineparams

def DbOpen():
	"""Opens the database and returns database connection."""
	table_creation_necessary=(not os.path.exists(pipelineparams.db_file_name))
	db_connection = sqlite3.connect(pipelineparams.db_file_name)
	if table_creation_necessary:
		# Create the run table.
		cur = db_connection.cursor()
		cur.execute("create table run (id primary key,run_directory,state)")
		db_connection.commit()
	return db_connection

def DbClose(db_connection):
	"""Closes the database."""
	db_connection.close()
