#!/usr/bin/env python
# create_samplesheet.py - given the name of a run folder, create the sample
# sheet for the run.

import os
import sys
import datetime as dt
import hcidemux.gnomex
import hcidemux.pipelineparams as params

def EraseCommas(s):
	"""Removes all commas in string s."""
	if s:
		return ''.join(s.split(','))
	else:
		return s

class MissingBarcodeError(LookupError):
      '''Barcode missing in lane where required for demultiplexing.'''

def CleanRow(rec, sample_count, dirname):
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
    lane = rec[1]
    if sample_count[lane] == 1:
        dr['index'] = ''
        dr['index2'] = ''
    elif sample_count[lane] > 1 and not len(dr['index']):
        raise MissingBarcodeError('Multiple samples in lane but missing barcode for sample: ' + dr['Sample_ID'])

    # Drop potential trailing integers on sequence request number 
    # e.g. 1480R1 -> 14806R
    # still works on clean sequence request number
    dr['Sample_Project'] = dr['Sample_Project'].split('R')[0] + 'R'

    row  = [dr['Lane'],dr['Sample_ID'],
            '_'.join([dr['Sample_ID'],dirname]),
            dr['Sample_Project'],
            dr['index'],dr['index2']]
    return(row)

class IEMSampleSheetCreator:
    def __init__( self, rundirectory, lanes=[] ):
        self.dirname = rundirectory
        self.id = self.dirname.split('/')[-1]
        self.lanes = lanes[:]
    def log( self, values ):
        sys.stderr.write(" ".join(values) + "\n")
    def CreateSampleSheet(self):
        """
        1. Query GNomEx DB. 
        2. Write out a config file (i.e. sample sheet) in the style of        
        Illumina Experiment Manager. 
        """
        # Select lanes from flow cell that are bar coded.
        g = hcidemux.gnomex.GNomExConnection()
        connection = g.GnConnect(params.db_user,params.db_password,asdict=False)
        c = connection.cursor()
        query = """select flowcell.barcode,
            flowcellchannel.number,
            sample.number,
            genomebuild.genomebuildname,
            sample.barcodesequence,
            appuser.firstname, 
            appuser.lastname,
            request.number,
            sample.barcodesequenceb
            from flowcell
            join flowcellchannel on flowcellchannel.idflowcell=flowcell.idflowcell
            join sequencelane on sequencelane.idflowcellchannel = 
            flowcellchannel.idflowcellchannel
            join sample on sequencelane.idsample = sample.idsample
            left outer join genomebuild on sequencelane.idgenomebuildalignto = 
            genomebuild.idgenomebuild
            join request on sequencelane.idrequest = request.idrequest
            join appuser on request.idappuser = appuser.idappuser
            where flowcellchannel.filename = '%s'\n""" % self.id
        if self.lanes:
          # Append lane selection to where clause.
          query += "and flowcellchannel.number in (%s)\n" % str(','.join(map(str,self.lanes)))
          query += "order by flowcellchannel.number, sample.number;"

        try:
          c.execute(query)
        except pymssql.OperationError:
			sys.stderr.write( query )
			raise
        results = c.fetchall()
		# Open file
		#ofs = open(samplesheet_fname,'w')
        ofs = sys.stdout
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
            row = CleanRow(rec, sample_count, self.dirname)
            ofs.write(','.join(row) + '\n')
         
        # Check if this is a rapid run. This will be the case
        # if the flow cell is registered as a single lane, but there
        # will be two lanes of data for it.
        if len(sample_count.keys())==1 and \
        int(sample_count.keys()[0]) == 1 and \
        os.path.exists(os.path.join(self.dirname,'Data','Intensities','L002')):    
            self.log(["Run", self.id, "is a rapid run. Duplicating rows for lane 2 in sample sheet",self.dirname,"."])
            # Write the results.
            nrecords = len(results) - 1
            for i,rec in enumerate(results):
                row = CleanRow(rec, sample_count)
                row[0] = 2
                ofs.write(','.join(row) + '\n')
                
        # Close database.
        connection.close()
		# Return file name.
        
if __name__ == "__main__":
	if len(sys.argv) < 2 or sys.argv[1] == '--help':
		sys.stderr.write("%s creates a sample sheet for an Illumina\n" % sys.argv[0] )
		sys.stderr.write("run folder from information in GNomEx and writes it to stdout.\n" )
		sys.stderr.write("Use: %s <run_folder> [lane# lane# ... ]\n" % sys.argv[0] )
		sys.exit(1)
	s = IEMSampleSheetCreator( sys.argv[1],map(int,sys.argv[2:]))
	s.CreateSampleSheet()
