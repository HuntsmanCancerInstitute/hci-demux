import os
import sys
import smtplib
#from email.MIMEMultipart import MIMEMultipart
#from email.MIMEText import MIMEText
#from email.MIMEImage import MIMEImage

from email.mime.application import MIMEApplication
from email.mime.multipart import MIMEMultipart
from email.mime.text import MIMEText

class Emailer:
	"""Base class for objects that need to email stuff."""

	def send_msg( self, from_addr, to_addr, subject, message, attachments=[], msgtype="plain" ):
		msg = MIMEMultipart()
		msg['Subject'] = subject
		msg['From'] = from_addr
		msg['To'] = ", ".join(to_addr)
		msg.attach(MIMEText(message))

		for f in attachments:
			fp=open(f,"rb")
			msg.attach(MIMEApplication(fp.read(),
				Content_Disposition='attachment; filename="%s"' % os.path.basename(f),
				Name=os.path.basename(f)))

		s = smtplib.SMTP('hci-mail.hci.utah.edu')

		# Send the message.
		retval = s.sendmail( from_addr, to_addr, msg.as_string() )
		s.quit()
		return retval

# Example of use:
def main():
	if len(sys.argv) != 4:
		usage()
		sys.exit(1)
	
	flowcell = sys.argv[1]
	attachments = sys.argv[2:]
	for file in attachments:
		if not os.path.exists(file):
			print("File '%s' not found." % file)
			sys.exit(1)
	
	from_addr = "sbsuser@hci-bio2.hci.utah.edu"
	to_addr = ['brett.milash@hci.utah.edu',
		'brian.dalley@hci.utah.edu',
	]
	subject = "Flow cell %s processing complete." % flowcell
	message_content = """Data processing for flow cell %s is complete.
The GERALD summary and sample tracking forms are attached.
If you no longer wish to receive this automated email message
please contact Brett Milash (brett.milash@hci.utah.edu).
""" % flowcell

	Emailer().send_msg( from_addr, to_addr, subject, message_content, attachments )
