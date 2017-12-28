"""
Fastq.py - classes for reading and writing fastq data to and from files.
"""
import sys
import subprocess

class Sequence:
	def __init__( self,name,sequence,quality):
		self.name = name
		self.sequence = sequence
		self.quality = quality
		self.spacer = "+"

class Reader:
	"""
	Reader objects read Fastq files and return Sequence objects. If the
	Fastq file is gzipped (and its name ends with .gz) the input will
	be gunzipped. The original file will not be changed.
	"""
	def __init__(self,filename):
		self.ifs=None
		self.p=None
		if filename.endswith(".gz"):
			child_args=["/bin/gunzip","-c",filename]
			child_dir="."
			self.p = subprocess.Popen(args=child_args,cwd=child_dir,stdout=subprocess.PIPE,stderr=subprocess.PIPE)
			self.ifs=self.p.stdout
		else:
			self.ifs=open(filename,"r")
	
	def close(self):
		#if self.p:
		#	self.p.terminate()
		self.ifs.close()
		self.ifs=None
	
	def __iter__(self):
		return self

	def next( self ):
		"""
		Returns the next read as a Sequence object.
		"""
		name = self.ifs.readline().strip()
		sequence = self.ifs.readline().strip()
		spacer = self.ifs.readline().strip()
		quality = self.ifs.readline().strip()
		if not ( name and sequence and spacer and quality ):
			raise StopIteration
		return Sequence(name,sequence,quality)
	
class Writer:
	"""
	Writer objects write Fastq data to a file. If the file name ends with
	".gz" the data will by gzipped on the fly.
	"""
	def __init__(self,filename):
		self.ofs=None
		self.p=None
		if filename.endswith(".gz"):
			child_args=["/bin/gzip"]
			child_dir="."
			self.p = subprocess.Popen(args=child_args,cwd=child_dir,stdin=subprocess.PIPE,stdout=open(filename,"w"))
			self.ofs=self.p.stdin
		else:
			self.ofs=open(filename,"w")
	
	def write(self,sequence_object):
		"""
		write writes a single Sequence class object to the Fastq
		file.
		"""
		self.ofs.write(sequence_object.name + "\n")
		self.ofs.write(sequence_object.sequence + "\n")
		self.ofs.write(sequence_object.spacer + "\n")
		self.ofs.write(sequence_object.quality + "\n")
	
	def close(self):
		"""
		close closes the Fastq file.
		"""
		self.ofs.close()

def test():
	"""
	Tests Reader and Writer classes by reading a gzipped fastq
	file and writing it to a gzipped fastq file.
	"""
	infile = sys.argv[1]
	outfile = sys.argv[2]
	ifs=Reader(infile)
	ofs=Writer(outfile)
	count = 0
	for sequence_object in ifs:
		ofs.write(sequence_object)
		count+=1
	ifs.close()
	ofs.close()
	print "Read %d sequences from %s and wrote them to %s." % ( count, infile, outfile )


if __name__ == "__main__":
	test()
