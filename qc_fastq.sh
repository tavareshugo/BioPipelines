#!/bin/bash

###########################
#
# QC fastq files
# performs basic quality 
# filtering and
# adapter trimming
#
# Hugo Feb 2016
#
###########################


#
# Set user arguments
#
# Usage message
usage()
{
cat << EOF

usage: $0 <options>

fastq quality control and filtering.

OPTIONS:
  -h      show this Help message.
  -o      Output directory.
  -p      file with list of Prefix fastq files names (see details).
  -1      suffix of read1 files.
  -2      suffix of read2 files.
  -f      Filtering options passed to cutadapt. Do not include input and 
          output files, nor Illumina adapter removal (which are 
          options passed by this script).
  -c      number of processing cores to use (equivalent to '-t' option of
          FastQC).


DETAILS:

This script processes several fastq files producing quality reports 
(using FastQC) from raw and filtered files. Quality filtering (using 
cutadapt) includes adapter removal of the Illumina TruSeq Universal 
and Indexed adapters). The script uses GNU parallel to parallelize the 
filtering across multiple cores, so please make sure you have it 
installed.

As an input the script takes a list of read file name's prefixes. 
For example, if there were paired-end reads from two samples, located in:
./read_directory/sample1_r1.fq.gz
./read_directory/sample1_r2.fq.gz
./read_directory/sample2_r1.fq.gz
./read_directory/sample2_r2.fq.gz

We could create a file called "fastq_prefix.txt", containing the following 
file name prefixes:
./read_directory/sample1
./read_directory/sample2

The suffix corresponding to read1 and read2 in this example would be 
'_r1.fq.gz' and '_r2.fq.gz', respectively.

Therefore, the command would be:
$0 -p fastq_list.txt -1 _r1.fq.gz -2 _r2.fq.gz -o ./output_dir/

The script will create two output directories named 'fastqc' and 
'filtered_reads'. The first contains the FastQC reports for raw and
filtered reads. The second contains the filtered reads, the log files 
from filtering software and two .csv files with information compiled 
from the log files (useful for plotting).

EOF
}

# Get options
while getopts "ho:p:1:2:f:c:" OPTION
do
  case $OPTION in
    h)  usage; exit 1;;
    o)  outdir=$OPTARG;;
    p)  file_list=$OPTARG;;
    1)  r1_suf=$OPTARG;;
    2)  r2_suf=$OPTARG;;
    f)  options=$OPTARG;;
    c)  threads=$OPTARG;;
    ?)  usage; exit;;
  esac
done

# Check that all options were passed
if [[ -z $outdir ]] || [[ -z $file_list ]] || [[ -z $r1_suf ]] || [[ -z $r2_suf ]] || [[ -z $options ]]
then
  printf "\n=========================\n ERROR: missing options\n=========================\n\n"
  usage
  exit 1
fi



#
# Create output directories
#
mkdir -p "$outdir/fastqc/raw"
mkdir -p "$outdir/fastqc/filtered"
mkdir -p "$outdir/filtered_reads"



#
# Prepare input and output file names
# 
in1=""  #input read1 files
in2=""  #input read2 files
out_log=""  #output cutadapt log files
out1=""  #output read1 files
out2=""  #output read2 files

while read f 
do
	in1="$in1 ${f}${r1_suf}"
	in2="$in2 ${f}${r2_suf}"
	out_log="$out_log ${outdir}/filtered_reads/${f##*/}.filtered.log"
	out1="$out1 ${outdir}/filtered_reads/${f##*/}.filtered.${r1_suf/#./}"
	out2="$out2 ${outdir}/filtered_reads/${f##*/}.filtered.${r2_suf/#./}"
done < "$file_list"



#
# FastQC on raw files
#
printf "Starting FastQC on raw files\n"
fastqc -t $threads -o $outdir/fastqc/raw/ $in1 $in2
printf "Finished FastQC on raw files\n"



#
# Cutadapt filtering
#
# Using GNU parallel
printf "Starting cutadapt with options:\n  cutadapt $options\n\n"

parallel --xapply \
cutadapt \
-a TruSeq_indexed=AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC \
-A TruSeq_universal_rc=AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATT \
$options \
-o {1} -p {2} \
{3} \
{4} '>' \
{5} \
::: $out1 ::: $out2 ::: $in1 ::: $in2 ::: $out_log

printf "Finished cutadapt"



#
# FastQC on filtered files
#
printf "Starting FastQC on filtered files\n"
fastqc -t $threads -o $outdir/fastqc/filtered $out1 $out2
printf "Finished FastQC on filtered files\n"


#
# Parse cutadapt output to single .csv file
#
# This is a convenient format for plotting and comparing samples
printf "sample,type,count\n" > $outdir/filtered_reads/cutadapt_read_stats.csv
printf "sample,type,count\n" > $outdir/filtered_reads/cutadapt_basepair_stats.csv

for f in $out_log
do
	name=$(basename $f)  # take the basename of the log file
	
	# Pipe to parse output to csv
	grep -A 5 "Total read pairs processed" $f | \
	sed 's/[:,]//g' | \
	sed 's/^  //g' | \
	sed 's/([^)]*)//g' | \
	sed 's/ \+ /,/g' | \
	sed 's/ $//g' | \
	sed "s/^/$name,/g" >> $outdir/filtered_reads/cutadapt_read_stats.csv
	
	grep "Total basepairs processed" $f | \
	sed 's/,//g' | \
	sed 's/ bp//g' | \
	sed 's/: \+/,/g' | \
	sed "s/^/$name,/" >> $outdir/filtered_reads/cutadapt_basepair_stats.csv
	
	grep "Total written" $f | \
	sed 's/,//g' | \
	sed 's/ bp//g' | \
	sed 's/([^)]*)//g' | \
	sed 's/ : \+/,/g' | \
	sed "s/^/$name,/" >> $outdir/filtered_reads/cutadapt_basepair_stats.csv
done

