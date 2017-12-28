#!/usr/bin/env python
# handle_email.py - script to provide email interface for PipelineAutomator.
#
# An entry in /etc/aliases directs email to stdin for this program.

import sys
import os
import time

#outdir=os.path.realpath(os.path.dirname(sys.argv[0]))
outdir="/tmp"
outfile=os.path.join(outdir,"email_message.%d" % int(time.time()))

ofs=open(outfile,'w')

#content=sys.stdin.read()
#while content:
#	ofs.write(content)
#	content=sys.stdin.read()
ofs.write(sys.stdin.read())
ofs.close()

