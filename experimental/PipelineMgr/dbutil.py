#!/usr/bin/env python3
import sys
try:
        import sqlite3
except ImportError:
        import sqlite as sqlite3

TAB="\t"

def Dump( conn ):
	c = conn.cursor()
	c.execute("select id, run_directory, state from run order by id")
	recs = c.fetchall()
	for rec in recs:
		print(TAB.join(list(rec)))

def AddRun( conn, run_id, run_directory ):
	cursor = conn.cursor()
	t = ( run_id, run_directory, 'new' )
	try:
		retval1 = cursor.execute("insert into run (id,run_directory,state) values (%s,%s)",t )
		retval2 = conn.commit()
	except sqlite3._sqlite.IntegrityError:
		print("Run %s is already in the database." % run_id)
		return
	print("retval1: %s, retval2: %s" % ( retval1, retval2))

def ChangeState( conn, run_id, new_state ):
	cursor = conn.cursor()
	t = ( new_state, run_id )
	retval1=cursor.execute("update run set state = %s where id = %s", t )
	retval2=conn.commit()
	print("retval1: %s, retval2: %s" % ( retval1, retval2))

def DeleteRun( conn, run_id ):
	cursor = conn.cursor()
	retval1=cursor.execute("delete from run where id = '%s'" % run_id )
	retval2=conn.commit()
	print("retval1: %s, retval2: %s" % ( retval1, retval2))

def usage():
	print("Use: %s --help | --state <run_id> <state_name> | --delete <run_id> | --add <run_id> <run_directory> | --dump" % sys.argv[0])

def main():
	# If no args or incorrect args, print the usage message.
	if len(sys.argv) < 2 or sys.argv[1] not in ['--help','--state','--delete','--add','--dump']:
		usage()
		sys.exit(1)
	
	if sys.argv[1] == '--help':
		usage()
		sys.exit(0)

	conn = sqlite3.connect("pipelinemgr.db")

	if sys.argv[1] == '--delete':
		DeleteRun( conn, sys.argv[2] )
	elif sys.argv[1] == '--add':
		AddRun( conn, sys.argv[2], sys.argv[3] )
	elif sys.argv[1] == '--dump':
		Dump( conn )
	elif sys.argv[1] == '--state':
		ChangeState( conn, sys.argv[2], sys.argv[3] )
	else:
		usage()

if __name__ == "__main__":
	main()
