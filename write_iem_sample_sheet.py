#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Jan 22 10:20:00 2018

@author: Chris Conley
"""



# import GNomEx query. 
import pickle
# write to file
import sys
# get the date
import datetime as dt



fh = open('/Users/bioinformatics/datastore/demux-dev/gnomex_paired_end_dual_index.obj', 'r')
#fh = open('/Users/bioinformatics/datastore/demux-dev/gnomex_paired_end_single_index.obj', 'r')
results = pickle.load(fh)
fh.close()


def EraseCommas(s):
	"""Removes all commas in string s."""
	if s:
		return ''.join(s.split(','))
	else:
		return s

def CleanRow(rec):
    # drop commas 
    dr = {'Lane' : EraseCommas(str(rec[1])),
                'Sample_ID' : EraseCommas(rec[2]),
                'Sample_Project' : EraseCommas(rec[7]),
                'index' : EraseCommas(rec[4]),
                'index2' : EraseCommas(rec[8]) }
    
    # replace nones with emptry string
    # no trailing white space 
    for k,v in dr.items(): 
        if v is None:
            dr[k] = ''
        dr[k] = dr[k].strip() 
            
    # If lane has a single sample, don't write its barcode.
    if sample_count[lane] == 1:
        dr['index'] = ''
        dr['index2'] = ''
    
    row  = [dr['Lane'],dr['Sample_ID'],
            dr['Sample_ID'],dr['Sample_Project'],
            dr['index'],dr['index2']]
    return(row)
    

# write to standard out
ofs = sys.stdout

# Write header.
header = ['FCID','Lane','SampleID','SampleRef','Index','Description','Control','Recipe','Operator','SampleProject']

raw_date = dt.date.today()
todays_date = '/'.join([str(raw_date.month),str(raw_date.day),str(raw_date.year)])

# Figure out instrument from RunInfo.xml later.
header = [('IEMFileVersion', '4'), 
          ('Date', todays_date),
          ('Workflow', 'GenerateFASTQ'),
          ('Application', 'HiSeq FASTQ Only'),
          ('Assay', 'TruSeq HT'),
          ('Description', ''),
          ('Chemistry', 'Amplicon')]
ofs.write('[Header]\n')
for k,v in header: 
    ofs.write(str(k) + ',' + str(v) + '\n')

ofs.write('\n[Settings]\n')
settings = [('Adapter',''),
            ('AdapterRead2','')]

ofs.write('\n[Data]\n')
data_cols = ['Lane','Sample_ID','Sample_Name','Sample_Project','index','index2']
ofs.write(','.join(data_cols) + '\n')

# Count how many samples are in each lane.
sample_count = {}
for rec in results:
	try:
		lane = rec[1]
	except KeyError:
		print rec
		raise
	try:
		sample_count[lane]+=1
	except KeyError:
		sample_count[lane] = 1

# Write the results.
nrecords = len(results) - 1
for i,rec in enumerate(results):
    row = CleanRow(rec)
    if i == nrecords:
        ofs.write(','.join(row))
    else:
        ofs.write(','.join(row) + '\n')
         
        
# Check if this is a rapid run. This will be the case
# if the flow cell is registered as a single lane, but there
# will be two lanes of data for it.
if len(sample_count.keys())==1 and \
 int(sample_count.keys()[0]) == 1 and \
 os.path.exists(os.path.join(self.dirname,'Data','Intensities','L002')):    
     self.log(["Run", self.id, "is a rapid run. Duplicating rows for lane 2 in sample sheet",samplesheet_fname,"."])
     row = CleanRow(rec)
     row[0] = 2
     ofs.write(','.join(row) + '\n')