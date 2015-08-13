#!/bin/bash

############################
# Indexing reference fasta #
# Hugo Dec 2014            #
############################

#
# Usage information
#
usage()
{
cat << EOF

This script indexes a reference fasta file for bowtie2, 
samtools, Picard tools and bwa.

usage: $0 fasta_prefix

  fasta_prefix    full path to the fasta file, excluding the file 
                  extension. The file must have .fa extension.

EOF
}

if [ $# -eq 0 ]
then
	usage
	exit 1
fi

#
# Get user input
#
prefix=$1 #the prefix of the reference fasta file (including path). Fasta file will be assumed to have .fa extension

#
# Run all indexers
#
bowtie2-build $prefix.fa $prefix &

samtools faidx $prefix.fa &

java -jar ~/bin/picard.jar CreateSequenceDictionary REFERENCE=$prefix.fa OUTPUT=$prefix.dict &

bwa index $prefix.fa &

wait
