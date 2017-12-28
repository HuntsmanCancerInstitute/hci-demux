#!/usr/bin/python
# EvaluateIndexing.py - generates barcode processing report.
# TODO:
#	Needs to handle flow cells with no barcode read.
#
#	Doesn't display sample name for non-barcoded samples and lane
#	status is incorrect.
#


import sys
import os
import re
import threading
import xml.sax.xmlreader
import xml.sax.handler
from logger import Logger

TAB = '\t'

class RunInfoHandler( xml.sax.handler.ContentHandler ):
	'''XML parser event handler to count number of bases sequenced
	in this run.'''
	def __init__( self ):
		self.numbases = 0

	def startElement(self,name,attrs):
		if name == 'Read' and attrs.getValue('IsIndexedRead') == 'N':
			self.numbases += int(attrs.getValue('NumCycles'))

class ReadCounter(threading.Thread,Logger):
	def __init__( self, command, semaphore, sample, lane, barcode, parent ):
		threading.Thread.__init__(self)
		self.command=command
		self.semaphore = semaphore
		self.sample = sample
		self.lane = lane
		self.barcode = barcode
		self.parent = parent
	
	def run( self ):
		self.semaphore.acquire()
		self.Log( '*** Starting '+ self.command)
		ifs = os.popen(self.command, 'r' )
		lines = ifs.readlines()
		ifs.close()
		linecount = int( ''.join(lines))
		numreads = linecount / 4
		self.parent.StoreReadCount( numreads, self.barcode, self.lane )
		self.semaphore.release()

def DisplayAsGb( bases ):
	'''Converts number of bases into number of gigabases (actually billions of bases).'''
	return '%.2f Gb' % ( bases / 1000000000.0 )

def DisplayAsMb( bases ):
	'''Converts number of bases into number of megabases (actually millions of bases).'''
	return '%.2f Mb' % ( bases / 1000000.0 )

def DisplayPercent( fraction, total ):
	try:
		return "%2.1f%%" % (100.0 * fraction / total)
	except ZeroDivisionError:
		return "0.0%"

def DisplayIntCommas( n ):
	s = `n`				# Number as a string.
	letters=list(s)			# Number as list of letters.
	letters.reverse()
	s_rev = ''.join(letters)
	substrs = re.findall('..?.?', s_rev ) # Groups of 3 digits.
	s_rev = ','.join(substrs)
	letters=list(s_rev)
	letters.reverse()
	retval = ''.join(letters)
	return retval

class IndexingEvaluator(Logger):
	def __init__( self, run_name, pipeline_version ):
		# Maximum number of simultaneous gunzip processes.
		self.maxthreads=4
		# Semaphore for controlling concurrent threads.
		self.sem = threading.Semaphore(self.maxthreads)
		# Thread objects for parallel execution.
		self.threads=[]
		self.pipeline_version = pipeline_version

		self.run_folder_name = os.path.realpath(run_name)
		if not os.path.exists(self.run_folder_name):
			self.Log("%s: Run folder %s not found." % ( sys.argv[0], self.run_folder_name ) )
			sys.exit(1)
		# Build the location of the demultiplexed folder.
		if self.pipeline_version == '1.7':
			subdir = os.path.join( "Data","Intensities","BaseCalls","Demultiplexed")
		else:
			subdir = "Unaligned"
		self.demultiplex_folder = os.path.join(self.run_folder_name,subdir)
		self.Log( "demultiplex_folder: %s" % self.demultiplex_folder )
		os.chdir(self.demultiplex_folder)
	
		self.readcount = {}	# Dictionary of read counts by lane, barcode.
		self.samplename = {}	# Dictionary of sample names by lane, barcode.
		self.numsamples = {1:0,2:0,3:0,4:0,5:0,6:0,7:0,8:0}	# Dictionary of number of samples for each lane.
		self.requester = {}	# who requested which sample.

	def CountBases( self ):
		'''Parses the RunInfo.xml file in the run directory to
		determine how many bases of sequence data produced per read.'''
		handler=RunInfoHandler()
		xml.sax.parse(os.path.join(self.run_folder_name,"RunInfo.xml"),handler)
		self.numbases = handler.numbases

	def StoreReadCount( self, numreads, barcode, lane ):
		'''Stores result of counting reads on a sample's file.'''
		try:
			self.readcount[lane][barcode] = numreads
		except KeyError:
			self.readcount[lane] = { barcode: numreads }

	def CreateReadCounter( self, sample_regexp, lane, barcode, requester="unknown" ):
		'''Creates and returns a thread object that will
		gunzip and count the lines in a file.'''
		# sample_date_sequencer_run_barcode_lane[_end].txt.gz
		fname_regexp = "(%s)_[0-9]*_[A-Z0-9]*_[0-9]*_[A-Z0-9-]*_%d(_1)?\.txt\.gz$" % (sample_regexp,lane)
		self.Log(fname_regexp)
		# Find file.
		fname = None
		for file in os.listdir("."):
			m = re.match(fname_regexp,file)
			if m:
				fname = file
				sample=m.groups()[0]
				self.requester[sample] = requester
				# Record the lane/barcode/sample combination.
				try:
					self.numsamples[lane] += 1
					self.samplename[lane][barcode] = sample
				except KeyError:
					self.numsamples[lane] = 1
					self.samplename[lane] = { barcode: sample }
				break
		if fname is None:
			return None

		self.Log(fname)
		# Create a thread to count the reads in the file.
		cmd = "gunzip -c %s | wc -l" % fname
		return ReadCounter( cmd, self.sem, sample, lane, barcode, self )

	def Setup( self, selected_lanes ):
		self.Log("EvaluateIndexing.Setup: setting up QC report for pipeline version %s, lanes %s." % ( self.pipeline_version, selected_lanes ) )
		if self.pipeline_version == '1.7':
			messages = self.Setup_17( selected_lanes )
		else:
			messages = self.Setup_18( selected_lanes )
		return messages

	def EstimateBadBarcodes( self, lane ):
		'''Returns an estimated number of bad barcode reads for
		a lane. Based on number of bad barcode files for the lane.
		4M * n - 2M.'''
		# Count the number of *.fastq.gz files in the bad barcode directory.
		bad_barcode_dir = 'Undetermined_indices/Sample_lane%d' % lane
		pattern="lane[0-9]_Undetermined_L00[0-9]_R1_[0-9]*.fastq.gz"
		filecount = 0
		if os.path.exists(bad_barcode_dir):
			for filename in os.listdir(bad_barcode_dir):
				if re.match(pattern,filename):
					filecount+=1
		# Each file contains 4 million bad reads. The last file is 
		# probably half full.
		if filecount == 0:
			return 0
		else:
			return filecount * 4000000 - 2000000

	def TooManyBadBarcodes( self, lane ):
		'''Determines if there are too many bad barcodes to count them all, or whether its better to
		estimate.'''
		# Count the number of *.fastq.gz files in the bad barcode directory.
		if self.EstimateBadBarcodes(lane) > 20000000:
			return True
		else:
			return False

	def Setup_18( self, selected_lanes):
		'''Setup creates threads to count the reads in each sample-specfic Fastq file, or to
		count the reads per barcode in the files containing all the unexpected barcodes (i.e.
		the <lane>_sequence.txt files in the unknown directory). These threads get run by
		the RunThreads method. Each thread is responsible for calling a method to store its own
		results.'''

		self.CountBases()
		messages = []

		# Process the created_samplesheet.csv file. This lists the known barcodes on the flow cell.
		# For each lane/barcode/sample, open the data file and count the records. Record the
		# number of reads for the lane/barcode, and record that the barcode is associated with
		# a sample. Also store the number of samples for the lane.
		ifs = open(os.path.join(self.run_folder_name,"Data/Intensities/BaseCalls/created_samplesheet.csv"))
		for rec in ifs:
			f = rec.strip().split(',')
			try:
				lane = int(f[1])
				if lane not in selected_lanes:
					self.Log(" Skipping lane %d, sample %s, barcode %s. Not in selected_lanes %s." % ( lane, f[2], f[4], selected_lanes ))
					continue
				sample = f[2].upper()
				barcode = f[4].strip()
				try:
					requester = f[5].strip()
				except IndexError:
					requester = 'unknown'
				self.Log("Lane %d, sample %s, barcode %s." % ( lane, sample, barcode ))
			except ValueError:
				continue
			# Locate the data file. Create a thread to count the reads in the file.
			t = self.CreateReadCounter( sample, lane, barcode, requester )
			if t is not None:
				self.threads.append( t )
			else:
				# Complain.

				message = "Problem! No data file found for sample %s, lane %d." % ( sample, lane )
				self.Log( message )
				messages.append( message )
		ifs.close()

		# Identify any lanes not covered in the SamplesDirectories.csv
		# file. These would be lanes with a single sample where barcode
		# processing wasn't necessary.
		s=set(selected_lanes)		# Set of all selected lanes.
		f=set(self.samplename.keys())	# Set of barcoded lanes.
		missing_lanes=s.difference(f)
		self.Log("*** Samples without bar codes in lane(s) "+`missing_lanes`)
		for lane in missing_lanes:
			# Create a read counter for this lane's file.
			sample = '[0-9]*X[0-9]*'
			barcode = 'None'
			t = self.CreateReadCounter( sample, lane, barcode )
			if t is not None:
				self.threads.append( t )
			else:
				# Complain.
				message = "Problem! No data file found for lane %d." % lane
				self.Log( message )
				messages.append( message )

		# Process the "Undetermined_indices" directory. The sequence files there contain the reads not 
		# associated with known barcodes for the lane.
		for lane in selected_lanes:
			self.Log("Counting bad barcode reads for lane %d." % lane)
			self.readcount[lane]={}
			# Determine if too many bad barcode reads to count. If so, estimate based on
			# a small subset of the bad barcodes.
			if self.TooManyBadBarcodes(lane):
				# Only process a few of the files.
				file_pattern = "Undetermined_indices/Sample_lane%d/lane%d_Undetermined_L00%d_R1_00[12345].fastq.gz" % (lane,lane,lane)
				message = "Too many bad bar codes to count in lane %d. Estimating." % lane
				increment = self.EstimateBadBarcodes(lane) / 20000000.0
				self.Log(message)
				messages.append( message )
			else:
				# Process all the files for the first read.
				file_pattern = "Undetermined_indices/Sample_lane%d/lane%d_Undetermined_L00%d_R1_*.fastq.gz" % (lane,lane,lane)
				increment = 1

			ifs = os.popen("gunzip -c " + file_pattern)
			for rec in ifs:
				if rec[0] != '@':
					continue
				if self.pipeline_version == '1.7':
					barcode = rec.strip().split('#')[1][0:6]
				else:
					barcode = rec.strip().split(':')[-1]
					if not barcode:
						barcode = 'none'

				try:
					try:
						self.readcount[lane][barcode]+=increment
					except KeyError:
						self.readcount[lane][barcode]=increment
				except MemoryError:
					message = "Barcode problem with lane %d - suspect incorrect barcodes listed for lane." % lane
					self.Log(message)
					messages.append( message )
					self.readcount[lane]={}
					break
			ifs.close()
		return messages

	def Setup_17( self, selected_lanes ):
		'''Setup creates threads to count the reads in each sample-specfic Fastq file, or to
		count the reads per barcode in the files containing all the unexpected barcodes (i.e.
		the <lane>_sequence.txt files in the unknown directory). These threads get run by
		the RunThreads method. Each thread is responsible for calling a method to store its own
		results.'''

		self.CountBases()
		messages = []

		# Process SamplesDirectories file. This lists the known barcodes on the flow cell.
		# For each lane/barcode/sample, open the data file and count the records. Record the
		# number of reads for the lane/barcode, and record that the barcode is associated with
		# a sample. Also store the number of samples for the lane.
		ifs = open("SamplesDirectories.csv")
		for rec in ifs:
			f = rec.strip().split(',')
			try:
				lane = int(f[1])
				if lane not in selected_lanes:
					continue
				sample = f[2].upper()
				barcode = f[4].strip()
				try:
					requester = f[5].strip()
				except IndexError:
					requester = 'unknown'
				self.Log("Lane %d, sample %s, barcode %s." % ( lane, sample, barcode ))
			except ValueError:
				continue
			# Locate the data file. Create a thread to count the reads in the file.
			t = self.CreateReadCounter( sample, lane, barcode, requester )
			if t is not None:
				self.threads.append( t )
			else:
				# Complain.
				self.Log( "Problem! No data file found for sample %s, lane %d." % ( sample, lane ) )
		ifs.close()

		# Identify any lanes not covered in the SamplesDirectories.csv
		# file. These would be lanes with a single sample where barcode
		# processing wasn't necessary.
		s=set(selected_lanes)		# Set of all selected lanes.
		f=set(self.samplename.keys())	# Set of barcoded lanes.
		missing_lanes=s.difference(f)
		self.Log("*** No barcoded samples in lane(s) "+`missing_lanes`)
		for lane in missing_lanes:
			# Create a read counter for this lane's file.
			sample = '[0-9]*X[0-9]*'
			barcode = 'None'
			t = self.CreateReadCounter( sample, lane, barcode )
			if t is not None:
				self.threads.append( t )
			else:
				# Complain.
				self.Log( "Problem! No data file found for lane %d." % lane )

		# Process the "unknown" directory. The sequence files there contain the reads not 
		# associated with known barcodes for the lane.
		# Find the most recent "GERALD*" directory in the unknown directory.
		pathname='unknown'
		l = os.listdir(pathname)
		l.sort()
		l.reverse()
		for subdir in l:
			if subdir[0:6] == 'GERALD':
				pathname = os.path.join(pathname,subdir)
				break
		# Locate the s_%d_sequence.txt.gz files. These are the fastq files with 
		# unexpected barcodes for each lane.
		for lane in selected_lanes:
			self.readcount[lane]={}
			for fname in [ "s_%d_sequence.txt.gz" % lane, "s_%d_1_sequence.txt.gz" % lane ]:
				if os.path.exists(os.path.join(pathname,fname)):
					fname = os.path.join(pathname,fname)
					self.Log( 'Processing ' + fname )
					ifs = os.popen("gunzip -c %s" % fname )
					for rec in ifs:
						if rec[0] != '@':
							continue
						barcode = rec.strip().split('#')[1][0:6]
						try:
							self.readcount[lane][barcode]+=1
						except KeyError:
							self.readcount[lane][barcode]=1
					ifs.close()
					break
		return messages

	def RunThreads( self ):
		# Start the threads to count the reads.
		for thread in self.threads:
			thread.start()
		# Wait for threads to finish.
		while self.threads:
			self.Log( "Waiting for %d of %d threads." % ( len(self.threads),len(self.threads)))
			self.threads[-1].join()
			self.threads.pop()

	def ReportHeader( self, ofs, messages=[] ):
		column_headers = [ 'Lane','Status','Sample','Index','% Reads','# Reads','Volume','Requester']
		ofs.write( 'Barcode Processing Summary:\t'+self.run_folder_name+'\n')
		ofs.write('\n')
		if messages:
			ofs.write('\n'.join(messages)+'\n')
		ofs.write( TAB.join(column_headers) + '\n' )

	def LaneSummary( self, lane, total_reads, total_bases, ofs ):
		'''Creates simple lane summary.'''
		ofs.write( TAB.join([ 'Lane %d summary'%lane,'','','','',DisplayIntCommas( total_reads ), DisplayAsGb( total_bases )]) + '\n' )

	def FlowcellSummary( self, total_reads, total_bases, ofs ):
		'''Creates simple flowcell summary.'''
		ofs.write( TAB.join([ 'Flowcell summary', '','','','',DisplayIntCommas( total_reads ), DisplayAsGb( total_bases )]) + '\n' )

	def GenerateReport( self, outfile="-", messages=[] ):
		"""Generates text report. Tab delim."""
		self.Log("Generating report.")
		if outfile=="-":
			ofs = sys.stdout
		else:
			ofs = open(outfile, 'w' )
		barcoded_lanes = self.readcount.keys()
		# Report header.
		self.ReportHeader( ofs, messages )

		total_reads=0
		total_bases=0
		# for each lane
		for lane in range(1,9):
			if lane in barcoded_lanes:
				(lane_reads,lane_bases)=self.GenerateLaneReportBarcoded( lane, ofs )
				self.LaneSummary(lane,lane_reads,lane_bases,ofs)
				total_reads+=lane_reads
				total_bases+=lane_bases
				ofs.write('\n')
			else:
				self.GenerateLaneReportNoBarcode( lane, ofs )

		# Flowcell Summary.
		self.FlowcellSummary( total_reads, total_bases, ofs )
		if ofs != sys.stdout:
			ofs.close()

	def GenerateLaneReportNoBarcode( self, lane, ofs ):
		pass

	def GenerateLaneReportBarcoded_Heap( self, lane, ofs ):
		self.Log("Generating report for lane %d (Heap implementation)." % lane )
		total_reads=0
		total_bases=0
		# Place holder for report rows before formatting.
		rows = []

		# Readcounts is a list of (barcode,count) tuples.
		readcounts = self.readcount[lane].items()
		# Create a list of (count,barcode,samplename) tuples.
		l = []
		lane_readcount_total = 0
		for ( barcode, count ) in readcounts:
			lane_readcount_total += int(count)
			try:
				l.append( ( int(count), barcode, self.samplename[lane][barcode] ) )
			except KeyError:
				l.append( ( int(count), barcode, None ) )

		# Create heap 

		# Sort by decreasing count.
		l.sort()
		l.reverse()

		# Determine if any unexpected barcodes found more frequently
		# than the known barcodes. This would indicate a problem!
		samples=map(lambda x: x[2], l[0:self.numsamples[lane]])
		note = 'OK'
		if None in samples:
			# Problem! An unexpected barcode present at higher frequency
			# than one or more expected barcodes.
			note = 'Problem!'

		# Determine if more than 10% of the barcodes are uninterpretable.
		# This would also be a problem.

		num_known_samples = 0
		i = 0
		while i < len(l) and num_known_samples < self.numsamples[lane]:
			(count,barcode,sample) = l[i]
			bases = count * self.numbases
			total_reads += count
			total_bases += bases
			percentage = DisplayPercent( count, lane_readcount_total )
			row = [ `lane`, note, sample or 'None', barcode, percentage, DisplayIntCommas(count), DisplayAsMb(bases), self.requester.get(sample,'unknown') ]
			ofs.write( TAB.join(row) + '\n' )
			if sample is not None:
				num_known_samples += 1
			i+=1

		# process rest of samples in lane.
		other_count = 0
		while i < len(l):
			(count,barcode,sample) = l[i]
			other_count += count
			i+=1
		try:
			other_ratio = other_count / lane_readcount_total
		except ZeroDivisionError:
			other_ratio = 0.0
		if other_ratio > 0.10:
			note = 'Problem!'
		other_bases = other_count * self.numbases
		percentage = DisplayPercent( other_count, lane_readcount_total )
		row = [ `lane`, note, 'None', 'others', percentage, DisplayIntCommas(other_count), DisplayAsMb(other_bases) ]
		ofs.write( TAB.join(row) + '\n' )
		total_reads+=other_count
		total_bases+=other_bases
		return (total_reads,total_bases)

	def GenerateLaneReportBarcoded( self, lane, ofs ):
			
		self.Log("Generating report for lane %d." % lane )
		total_reads=0
		total_bases=0
		# Place holder for report rows before formatting.
		rows = []

		# Readcounts is a list of (barcode,count) tuples.
		readcounts = self.readcount[lane].items()
		# Create a list of (count,barcode,samplename) tuples.
		l = []
		lane_readcount_total = 0
		for ( barcode, count ) in readcounts:
			# Casting as an int because count may not be an integer
			# due to estimating # of bad barcode reads.
			lane_readcount_total += int(count)
			try:
				l.append( ( int(count), barcode, self.samplename[lane][barcode] ) )
			except KeyError:
				l.append( ( int(count), barcode, None ) )
		# Sort by decreasing count.
		l.sort()
		l.reverse()

		# Determine if any unexpected barcodes found more frequently
		# than the known barcodes. This would indicate a problem!
		samples=map(lambda x: x[2], l[0:self.numsamples[lane]])
		note = 'OK'
		if None in samples:
			# Problem! An unexpected barcode present at higher frequency
			# than one or more expected barcodes.
			note = 'Problem!'

		# Determine if more than 10% of the barcodes are uninterpretable.
		# This would also be a problem.

		samplenames = self.samplename[lane].values()
		num_known_samples = 0
		i = 0
		while i < len(l) and num_known_samples < self.numsamples[lane]:
			(count,barcode,sample) = l[i]
			bases = count * self.numbases
			# Break out of loop if we're down to a fraction of a
			# percent of the lane's reads.
			if count < lane_readcount_total * 0.001:
				break
			total_reads += count
			total_bases += bases
			percentage = DisplayPercent( count, lane_readcount_total )
			row = [ `lane`, note, sample or 'None', barcode, percentage, DisplayIntCommas(count), DisplayAsMb(bases), self.requester.get(sample,'unknown') ]
			ofs.write( TAB.join(row) + '\n' )
			if sample is not None:
				num_known_samples += 1
				samplenames.remove(sample)
			i+=1

		# Report on any samples that haven't been documented yet.
		# This is necessary because samples with very few reads (<1%
		# of total) won't have been listed yet.
		if samplenames:
			for j in range(i,len(l)):
				(count,barcode,sample) = l[j]
				# If we found an actual sample...
				if sample in samplenames:
					percentage = DisplayPercent( count, lane_readcount_total )
					bases = count * self.numbases
					row = [ `lane`, note, sample or 'None', barcode, percentage, DisplayIntCommas(count), DisplayAsMb(bases), self.requester.get(sample,'unknown') ]
					ofs.write( TAB.join(row) + '\n' )
					samplenames.remove(sample)

		# If any samples remain in the list at this point they have
		# 0 reads.
		for sample in samplenames:
			# Find its barcode. Doesn't need to be efficient - this
			# should rarely happen.
			for barcode in self.samplename[lane].keys():
				if self.samplename[lane][barcode] == sample:
					break
			count = 0
			percentage = DisplayPercent( count, lane_readcount_total )
			bases = count * self.numbases
			row = [ `lane`, note, sample or 'None', barcode, percentage, DisplayIntCommas(count), DisplayAsMb(bases), self.requester.get(sample,'unknown') ]
			ofs.write( TAB.join(row) + '\n' )

		# process rest of samples in lane.
		other_count = 0
		while i < len(l):
			(count,barcode,sample) = l[i]
			other_count += count
			i+=1
		try:
			other_ratio = other_count / lane_readcount_total
		except ZeroDivisionError:
			other_ratio = 0.0
		if other_ratio > 0.10:
			note = 'Problem!'
		other_bases = other_count * self.numbases
		percentage = DisplayPercent( other_count, lane_readcount_total )
		row = [ `lane`, note, 'None', 'others', percentage, DisplayIntCommas(other_count), DisplayAsMb(other_bases) ]
		ofs.write( TAB.join(row) + '\n' )
		total_reads+=other_count
		total_bases+=other_bases
		return (total_reads,total_bases)

def usage():
	sys.stderr.write("Usage: %s <run_folder_name> <output_file_name> <pipeline_version>\nWhere pipeline version is '1.7' or '1.8'\n" % sys.argv[0] )

def RunReport( runname, outputfile, pipeline_version, selected_lanes=[1,2,3,4,5,6,7,8] ):
	e = IndexingEvaluator( runname, pipeline_version )
	messages = e.Setup(selected_lanes)
	e.RunThreads()
	e.GenerateReport( outputfile, messages )

def main():
	if len(sys.argv) != 4:
		usage()
		sys.exit(1)
	RunReport( sys.argv[1], sys.argv[2], sys.argv[3] )

if __name__ == "__main__":
	main()
