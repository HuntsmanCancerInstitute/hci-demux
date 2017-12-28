# hci-demux

hci-demux is a software package that automatically configures demultiplexing settings,
calls the Illumina BCL to Fasta conversion software, distributes data files to GNomEx,
and generates demultiplexing QC reports.

## Installation

To install the software:
1. clone the repository somewhere under the home directory of the account that will
run the pipelines.
2. Make sure bcl2fastq version 1.8.4 is installed.
3. Load the python/2.7.3 module (or any python 2.7 module) in the environment of the account running the software.
4. Make sure that python has installed pymssql and sysv_ipc libraries.
5. Edit the process_pipelines.sh script at the root of the repository if necessary to
redirect log files somewhere besides $HOME/Pipeline/logs.
6. Check the pipelineparams.py file out of Confluence and edit it appropriately. Place this file into
the hcidemux directory underneath the root of the project.
7. Edit the crontab script to reflect the install location of the software and the time at which you'd like
the software to run. Add the lines in the pipeline.crontab file to the crontab of the account running the
pipelines.

## Logging

The hourly crontab output will be sent to the software root directory, and this file will get overwritten each hour.
Detailed daily log files will be created in $HOME/Pipeline/logs, and these files are not overwritten. You may want to 
add a crontab entry to remove logfiles older than some number of days.
