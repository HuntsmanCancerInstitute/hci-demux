This directory contains code for the next iteration of the pipeline
automation. This will provide for a variety of demultiplexing options
which are now handled by shell scripts. 

Demultiplexing options will include:
	single-base mismatch using all barcode reads (current default option)
	molecular barcode processing
	no demultiplexing
	no mismatch demultiplexing

Run types:
	single end
	paired end

Bar code read options:
	single bar code read
	dual bar code read
	no bar code read

Sequencer options:
	HiSeq
	MiSeq

This system will be run by cron, but will also allow manual runs, and will
have a limited email interface to permit re-run triggering by the
sequencing lab staff.

Brett Milash, 05/09/2017.

Top-level program - executed via cmd line or crontab:
	exclusivity via semaphore.
	create db if needed.
	delete nonexistent runs from db.
	add new runs to db.
	Process existing runs. For each run:
		execute code attached to currrent state, update state
			as needed until no state change.
	States:
		new 
		waiting to be registered
		waiting for data
		pipeline
		reprocess
		complete
		error
	Simplified pipeline state - active execution of sample sheet
		creation, base calling, demultiplexing, file assembly,
		distribution to GNomEx, QC report generation, cleanup,
		(and reprocessing if needed). Each flow cell processed
		by one or more RunProcessors, which encapsulate complexity
		of number of reads, bar code reads, bar code lengths, 
		demultiplexing details, etc. A RunProcessor handles one
		or more lanes with identical characteristics
	Run class - a Run object represents the results from one run of 
		the sequencer, held in one run directory. A Run is processed
		by one or more RunProcessors.
	

Db - stores runs and current state. Single table:
	run id
	run directory
	state
