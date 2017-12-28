import time
import sys

class Logger(object):
	"""Base class for objects that need to log stuff."""
	def Log( self, msg ):
		lt = time.localtime( time.time() )
		timestamp = time.strftime("%Y-%m-%d %H:%M ",lt)
		if type(msg) == type([]):
			sys.stderr.write(timestamp+" ".join(map(str,msg)) + "\n")
		else:
			sys.stderr.write(timestamp+str(msg)+"\n")
		sys.stderr.flush()
