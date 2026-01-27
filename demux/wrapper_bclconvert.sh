#!/bin/bash

shopt -s extglob
set -euo pipefail

HEADER='
\n# Script Name: wrapper_bclconvert.sh
\n# Description: Wrapper script to demux NovaSeq 6000/X runs using bcl-convert.
\n# Version: 1.3.3 (2026-01-20)
\n# Author: Filipe G. Vieira
\n# Mail: fgvieira@sund.ku.dk
'
BASEDIR=`dirname $0`


module load python/3.12.8
# Test python imports
python3 -c 'import samshee'
python3 -c 'import argparse'
python3 -c 'import logging'
python3 -c 'import numpy'
python3 -c 'import pandas'
python3 -c 'import pathlib'
python3 -c 'import collections'
python3 -c 'import plotly'
python3 -c 'import matplotlib'


## Functions
function check_se() {
    RET=`zcat $1 | awk --field-separator "\t" 'NR%4==1 && !/^@/' | head --lines 1`
    if [[ $RET ]]; then
	echo -e "ERROR: malformed SE FASTQ file on line:\n$RET"
	exit 2
    fi
}

function check_pe() {
    RET=`paste <(zcat $1 | cut --delimiter " " --fields 1) <(zcat $2 | cut --delimiter " " --fields 1) | awk --field-separator "\t" 'NR%4==1 && ($1!=$2 || !/^@/)' | head --lines 1`
    if [[ $RET ]]; then
	echo -e "ERROR: malformed PE FASTQ file on line:\n$RET"
	exit 2
    fi
}

function ss_pool() {
    python3 -c "from samshee.samplesheetv2 import read_samplesheetv2; ss = read_samplesheetv2('$1'); [print(val, tag.removeprefix('PoolLane')) for tag, val in ss.header.items() if tag.startswith('PoolLane')]" | datamash -t " " --output-delimiter=":" groupby 1 collapse 2
}

function ss_proj() {
    python3 -c "from samshee.samplesheetv2 import read_samplesheetv2; ss = read_samplesheetv2('$1'); [print(data['Sample_Project']) for data in ss.applications['BCLConvert']['data']]" | sort -u
}
export -f check_se check_pe ss_pool ss_proj


## Setup
THREADS=5

IN_FOLDER=$1; shift
SS=$1; shift
OUT_FOLDER=$1; shift
EXTRA=$@

RUN=20`basename $IN_FOLDER`
if [[ `realpath $OUT_FOLDER` == "/maps/datasets/caeg_fastq" ]]; then
    OUT_FOLDER=$OUT_FOLDER/${RUN:0:4}/$RUN
else
    OUT_FOLDER=$OUT_FOLDER/$RUN
fi

# Check if output folder exists
if [ -d $OUT_FOLDER ]; then
    RND_STR=`mktemp --directory --dry-run`
    RND_STR=`basename $RND_STR`
    OUT_FOLDER=$OUT_FOLDER-$RND_STR
fi
mkdir -p $OUT_FOLDER


{
    # Log script info
    echo -e $HEADER

    ## Demultiplex
    echo `date`" [$RUN] Demultiplexing from $IN_FOLDER to $OUT_FOLDER with SampleSheet $SS (and extra: $EXTRA)"
    bcl-convert --bcl-input-directory $IN_FOLDER --output-directory $OUT_FOLDER --sample-sheet $SS --bcl-sampleproject-subdirectories true --force $EXTRA

    ## Get projects from SS
    PROJS=(`ss_proj $SS`)
    ## Check FASTQ files
    cd $OUT_FOLDER/
    for PROJ in ${PROJS[*]}
    do
	cd $PROJ/
	## Run gzip integrity in parallel
	echo `date`" [$RUN][$PROJ] Checking GZip files integrity"
	ls *.fastq.gz | xargs --max-procs $THREADS --delimiter '\n' --max-args 1 gzip --test


	## Check if FASTQ is valid
	#echo `date`" [$RUN][$PROJ] Checking if FASTQ files are valid"
	FASTQ_R1=(*_R1_001.fastq.gz)
	FASTQ_R2=(*_R2_001.fastq.gz)
	# https://unix.stackexchange.com/questions/651119/parallelize-a-function-using-xargs-and-separate-variables
	if [[ ${#FASTQ_R1[@]} -eq ${#FASTQ_R2[@]} ]]; then
	    echo `date`" [$RUN][$PROJ] Paired-End run"
	    #ls *.fastq.gz | xargs --max-procs $THREADS --delimiter '\n' --max-args 2 bash -c 'check_pe $1 $2' bash
	elif [[  ${#FASTQ_R2[@]} -eq 1 ]]; then
	    echo `date`" [$RUN][$PROJ] Single-End run"
	    #ls *.fastq.gz | xargs --max-procs $THREADS --delimiter '\n' --max-args 1 bash -c 'check_se $1' bash
	else
	    echo `date`" [$RUN][$PROJ] Error inferring if PE or SE sequencing"
	    exit 1
	fi


	## Create md5sum file
	echo `date`" [$RUN][$PROJ] Creating MD5 checksums"
	ls *.fastq.gz | xargs --max-procs $THREADS --delimiter '\n' --max-args 1 md5sum > $RUN.$PROJ.md5


	# Check duplicated MD5 checksums
	echo `date`" [$RUN][$PROJ] Checking for duplicated md5 checksums"
	cut --delimiter " " --fields 1 $RUN.$PROJ.md5 | sort | uniq -d


	# Check FASTQ file sizes
	echo `date`" [$RUN][$PROJ] Checking FASTQ files with size < 1Mb"
	find . -name "*.fastq.gz" -size 1M -printf "%k Kb\t%p\n" | grep -Fv Undetermined || true


	# Exit PROJ
	cd ../
    done

    ## Get pools from SS
    POOLS=(`ss_pool $SS`)
    ## Check cross-contamination
    for POOL in ${POOLS[*]}
    do
	echo `date`" [$RUN][$POOL] Check cross-contamination"
	python3 $BASEDIR/cross_contamination.py --index-counts $OUT_FOLDER/Reports/Index_Hopping_Counts.csv --index-known $BASEDIR/eDNA_index_list_UDP097-UDP288_UDI001-UDI096_250807.txt --lanes ${POOL#*:} --rpm-warn 100 --out-prefix $OUT_FOLDER/Reports/Index_Hopping_Counts/${POOL%:*}
    done

    TIMESTAMP=`date "+%Y%m%d_%H%M%S"`
    touch seqcenter.$TIMESTAMP.done
    cd ../
} 2>&1 | tee $OUT_FOLDER/$RUN.demux.log

exit 0
