#!/bin/bash
# process_pipelines.sh - part of hci-demux package to automate sequence
# processing.

# Create log directory if needed.
mkdir -p $HOME/Pipeline/logs

# Execute the demultiplexing pipelines. This will create the run database
# if it doesn't exist already.
python ./execute_pipeline.py >> `date +$HOME/Pipeline/logs/pipeline.log.\%Y-\%m-\%d` 2>&1

# Clean out the database of any runs that don't exist anymore.
for dir in `./dbutil.py --dump | cut -f2`
do
	#echo $dir
	runname=`basename $dir`
	#echo $runname
	if [ ! -d $dir ]
	then
		echo "Run directory $dir doesn't exist. Deleting it from database."
		./dbutil.py --delete $runname
	fi
done
