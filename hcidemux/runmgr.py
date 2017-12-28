# Class that manages all the runs in a specified .db file
import sys
import os
import re
import traceback
import sqlite3

from run import Run
from states import States
from logger import Logger
import pipelineparams as params
        
class RunMgr(Logger):
	def __init__(self):
		self.db_file = params.db_file
		
		# List of run objects that will get processed.
		self.runs_to_process = []
		
		# Open pipelinemgr database that tracks run processing. 
		# If it doesn't exist, a new database is created
		if self.DbExists():
			self.DbOpen()
		else:
			self.DbCreate()
		
	def AddRun( self, run_object ):
		self.runs_to_process.append( run_object )

	def DbExists(self):
			#Checks is piplinemgr database file exists
		return os.path.exists(self.db_file)
		
	def DbCreate(self):
		#create and open database
		self.Log("Creating database.")
		self.DbOpen()
		c = self.db_connection.cursor()
		
		# Creates runs table.
		c.execute("create table run (id primary key,dirname,state)")
		self.db_connection.commit()
		
		c.close()
		self.Log("Database creation complete.")

	def DbOpen( self ):
		"Open the database."""
		#self.Log("Opening database.")
		self.db_connection = sqlite3.connect(self.db_file)
		# Set up the database connection so the returned values
		# are byte strings rather than unicode.
		self.db_connection.text_factory = bytes
		
	def CleanUp( self ):
		"""Clean any old and deleted runs out of the database. Any
		Run in the database that doesn't appear on the filesystem
		can be removed from the database."""
		# Get all the runs out of the database.
		runs = self.GetAllRuns()
		c = self.db_connection.cursor()
		# For each run...
		for (id,dirname,state) in runs:
			# ... if not present on the file system...
			if not os.path.exists(dirname):
				# ... remove the record from the database.
				self.Log("Deleting old run %s (%s) from database." % (id,dirname))
				c.execute("delete from run where id = ?",(id,))
				self.db_connection.commit()

	def GetActiveRuns( self ):
		"""Returns list of active runs from database."""
		c = self.db_connection.cursor()
		c.execute("select id, dirname, state from run where state != ? and state != ?",(States.error_notified,States.complete))
		recs = c.fetchall()
		return recs
	
	def GetAllRuns( self ):
		"""Return list of all known runs."""
		c = self.db_connection.cursor()
		c.execute("select id,dirname,state from run")
		recs = c.fetchall()
		return recs
	
	def GetAllRunIds( self ):
		"""Return list of all known run ids from database."""
		c = self.db_connection.cursor()
		c.execute("select id from run")
		recs = c.fetchall()
		runs = map( lambda x:x[0], recs )
		return runs

	def Discover( self, root_directories ):
		"""Look for any new runs. Enter them into the database
		in the starting state."""
		self.Log("Discovering new sequencing runs.")
		# Find subdirectories in the runs directory whose name matches the
		# regular expression.
		known_run_ids = self.GetAllRunIds()
		#self.Log(["Known run ids",known_run_ids])

		new_runs = []
		# Pattern consists of Date_SequencerName_RunNumber_Flowcell.
		# Pattern updated to match both HiSeq and MiSeq runs.
		# Brett Milash 3/20/2015.

		# Use this pattern to find both Hiseq and Miseq runs:
		pattern = re.compile("[0-9]*_[A-Z0-9]*_[0-9]*_[A-Z0-9-]*$")

		# Use this pattern to find only Hiseq runs:
		#pattern = re.compile("[0-9]*_[A-Z0-9]*_[0-9]*_[AB][A-Z0-9]*")

		for runs_dir in root_directories:
			for subdir in os.listdir(runs_dir):
				run_full_path=os.path.join(runs_dir,subdir)
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
						self.Log(["Subdirectory",subdir,"already known."])
		if new_runs:
			for (subdir,run_full_path) in new_runs:
				run = Run( subdir,run_full_path )
				run.DbAdd( self.db_connection )
				self.Log(["Discovered run",subdir])
		else:
			self.Log("No new sequencing runs.")
			
	def ProcessRuns(self):
		# runs_to_process is set different for each pipeline called
		# set by giving RunMgr.runs_to_process = list
		# list is list of runs
		self.runs_to_process.sort()
		for run in self.runs_to_process:
			try:
				prev_state = None
				while run.state!= prev_state:
					prev_state = run.state
					try:
						if not run.state in run.transition:
							self.Log(["Skipping run",run.id,", no transitions from state", run.state])
							break
						(method,success_state,fail_state) = run.transition[run.state]
					
						if method(run):
							# Function successful
							self.Log(["Run",run.id,"transitioned from",run.state,"to",success_state])
							# Move run to new State.
							run.state=success_state
						else:
							self.Log(["Run",run.id,"transitioned from",run.state,"to",fail_state])
							run.state = fail_state
						run.DbUpdate(self.db_connection)
					except:
						(etype,evalue,etraceback) = sys.exc_info()
						details = traceback.format_exception( etype,evalue,etraceback)
						self.Log(["Exception encountered processing",run.id,". Continuing with next run. Details:" ]+details)
			except:
				run.state='o_error_detected'
				(etype,evalue,etraceback) = sys.exc_info()
				details = traceback.format_exception( etype,evalue,etraceback)
				self.Log(["Exception encountered processing",run.id,". Continuing with next run. Details:" ]+details)
				if run.NotifyError():
					run.state='p_error_notified'
				run.DbUpdate(self.db_connection)
		self.CleanUp()
