#!/bin/bash

###########
# Mapping pipeline using bowtie2
# Hugo Jul 2015
###########

# For debugging:
# bash mapping_pipeline.sh -o testoutdir -p testprefix -1 read1 -2 read2 -i testid -l testlibrary -b testbarcode -s testsample -r testref -c 6

#
# Picard and GATK paths - change the paths if necessary
#
picard="java -jar ~/bin/picard.jar"
gatk="java -jar ~/bin/GenomeAnalysisTK.jar"

#
# Set user arguments
#
# Usage message
usage()
{
cat << EOF
usage: $0 <options>

This script runs a mapping pipeline using bowtie2 (mapping), 
Picard tools (reorder and markdup) and GATK (realign around indels).
It also outputs several basic statistics for the final alignment.

OPTIONS:
  -h      Show this message
  -o      output directory
  -p      prefix for output files
  -1      path to read 1 fastq file
  -2      path to read 2 fastq file
  -i      sample ID (usually take the ID from fastq file header)
  -l      library name
  -b      library barcode (usually can be obtained from fastq header)
  -s      custom sample name
  -r      path to indexed reference genome (without file extension)
  -c      number of processing cores to use      

EOF
}

# Get options
while getopts “ho:p:1:2:i:l:b:s:r:c:” OPTION
do
  case $OPTION in
    h)  usage; exit 1;;
    o)  outdir=$OPTARG;;
    p)  outprefix=$OPTARG;;
    1)  read1=$OPTARG;;
    2)  read2=$OPTARG;;
    i)  id=$OPTARG;;
    l)  lb=$OPTARG;;
    b)  pu=$OPTARG;;
    s)  sm=$OPTARG;;
    r)  ref=$OPTARG;;
    c)  threads=$OPTARG;;
    ?)  usage; exit;;
  esac
done

# Check that all options were passed
if [[ -z $outdir ]] || [[ -z $outprefix ]] || [[ -z $read1 ]] || [[ -z $id ]] || [[ -z $lb ]] || [[ -z $pu ]] || [[ -z $sm ]] || [[ -z $ref ]] || [[ -z $threads ]]
then
  printf "\n=========================\n ERROR: missing options\n=========================\n\n"
  usage
  exit 1
fi


#
# Make output directory
#
mkdir -p $outdir/stats



#
# mapping with bowtie2 and coordinate sorting
#
bowtie2 -N 1 -p $threads \
-x $ref \
-1 $read1 -2 $read2 \
--rg-id $id \
--rg SM:$sm --rg PL:illumina --rg LB:$lb --rg PU:$pu \
|
samtools sort -T $outdir/samtools_tempfiles \
-o $outdir/temp1_$outprefix.bwt2.bam -;


#
# reorder reads to have same order as the reference
#
$picard ReorderSam \
INPUT=$outdir/temp1_$outprefix.bwt2.bam \
OUTPUT=$outdir/temp2_$outprefix.bwt2.reorder.bam \
REFERENCE=$ref.fa \
QUIET=true \
VERBOSITY=ERROR;

#rm $outdir/temp1*;


#
# remove duplicates with Picard MarkDuplicates
#
$picard MarkDuplicates \
INPUT=$outdir/temp2_$outprefix.bwt2.reorder.bam \
OUTPUT=$outdir/temp3_$outprefix.bwt2.reorder.markdup.bam \
METRICS_FILE=$outdir/stats/$outprefix.bwt2.reorder.markdup_metrics \
MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=1000 \
CREATE_INDEX=true

#rm $outdir/temp2*;


#
# Realign around indels
#
# Create target realigner
$gatk -T RealignerTargetCreator \
-nt $threads -R $ref.fa \
-o $outdir/temp3_$outprefix.bwt2.reorder.markdup.intervals \
-I $outdir/temp3_$outprefix.bwt2.reorder.markdup.bam

# Realign
$gatk -T IndelRealigner \
-R $ref.fa \
--baq CALCULATE_AS_NECESSARY \
-I $outdir/temp3_$outprefix.bwt2.reorder.markdup.bam \
-o $outdir/$outprefix.bwt2.reorder.markdup.rlgn.bam \
-targetIntervals $outdir/temp3_$outprefix.bwt2.reorder.markdup.intervals \
-noTags

#rm $outdir/temp3*


##### Collect statistics #####
#
# collect insert size distribution
#
$picard CollectInsertSizeMetrics \
INPUT=$outdir/$outprefix.bwt2.reorder.markdup.rlgn.bam \
OUTPUT=$outdir/stats/$outprefix.bwt2.reorder.markdup.rlgn.insert_size.hist \
HISTOGRAM_FILE=$outdir/stats/$outprefix.bwt2.reorder.markdup.rlgn.insert_size.pdf \
REFERENCE_SEQUENCE=$ref.fa &


#
# samtools basic stats
#
samtools flagstat $outdir/$outprefix.bwt2.reorder.markdup.rlgn.bam > $outdir/stats/$outprefix.bwt2.reorder.markdup.rlgn.flagstat &


#
# gatk callable sites stats
#
$gatk -T CallableLoci \
-I $outdir/$outprefix.bwt2.reorder.markdup.rlgn.bam \
-R $ref.fa \
--summary $outdir/stats/$outprefix.bwt2.reorder.markdup.rlgn.callable.summary \
-o $outdir/stats/$outprefix.bwt2.reorder.markdup.rlgn.callable.bed \
--format BED \
--minDepth 10 &

wait


