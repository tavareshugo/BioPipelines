#!/bin/bash

###########
# Mapping pipeline script
# Hugo Jul 2015
###########

# For debugging:
# bash mapping_pipeline.sh -o testoutdir -p testprefix -1 read1 -2 read2 -i testid -l testlibrary -b testbarcode -s testsample -r testref -c 6

#
# Picard and GATK paths - change the paths if necessary
#
picard="$HOME/bin/picard.jar"
gatk="$HOME/bin/gatk.jar"

#
# Set user arguments
#
# Usage message
usage()
{
cat << EOF
usage: $0 <options>

This script runs a mapping pipeline using either bowtie2 or bwa (mapping), 
Picard tools (reorder and markdup) and GATK (realign around indels).
It also outputs several basic statistics for the final alignment.

OPTIONS:
  -h      Show this message.
  -m      Mapping software to use. Choose between 'bwa' (default) or 'bowtie2'.
  -o      output directory.
  -p      prefix for output files.
  -1      path to read 1 fastq file.
  -2      path to read 2 fastq file.
  -i      sample ID (usually take the ID from fastq file header).
  -l      library name.
  -b      library barcode (usually can be obtained from fastq header).
  -s      custom sample name.
  -r      path to indexed reference genome (including file extension).
  -x      extra options passed to the mapper. This excludes read inputs,
          reference fasta, number of threads to use and read group info.
  -c      number of processing cores to use.
  -d      mark duplicates. Choose between 'yes' (default) or 'no'.
  -k      keep intermediate .bam files (these have a temp*_ prefix). 
          Choose between 'yes' or 'no' (default).

EOF
}

# Get options
while getopts “hm:o:p:1:2:i:l:b:s:r:x:c:d:k:” OPTION
do
  case $OPTION in
    h)  usage; exit 1;;
    m)  mapper=$OPTARG;;
    o)  outdir=$OPTARG;;
    p)  outprefix=$OPTARG;;
    1)  read1=$OPTARG;;
    2)  read2=$OPTARG;;
    i)  id=$OPTARG;;
    l)  lb=$OPTARG;;
    b)  pu=$OPTARG;;
    s)  sm=$OPTARG;;
    r)  ref=$OPTARG;;
    x)  extraopts=$OPTARG;;
    c)  threads=$OPTARG;;
    d)  dups=$OPTARG;;
    k)  keep=$OPTARG;;
    ?)  usage; exit;;
  esac
done


#
# Check options
#
# Check that all required options were passed
if [[ -z $outdir ]] || [[ -z $outprefix ]] || [[ -z $read1 ]] || [[ -z $id ]] || [[ -z $lb ]] || [[ -z $pu ]] || [[ -z $sm ]] || [[ -z $ref ]] || [[ -z $threads ]]
then
  printf "\n=========================\n ERROR: missing options\n=========================\n\n"
  usage
  exit 1
fi

# Define bowtie2 as the default mapper
if [[ $mapper = "" ]]
then
	mapper="bwa"
fi

if [[ $mapper != "bowtie2" ]] && [[ $mapper != "bwa" ]]
then
	printf "\nMapper has to be 'bowtie2' or 'bwa'\n" 1>&2
	exit 1
fi

# Define default removal of temp files
if [[ $keep = "" ]]
then
	keep="no"
fi

# Define default marking of duplicates
if [[ $dups = "" ]]
then
	dups="yes"
fi


#
# Make output directory
#
echo "Mapping pipeline outputs are in: $outdir" 1>&2
mkdir -p $outdir/stats


#
# Define filename of output
#
if [[ $dups == "yes" ]]
then
	output="$outprefix.$mapper.sort.reorder.markdup.rlgn"
elif [[ $dups == "no" ]]
then
	output="$outprefix.$mapper.sort.reorder.rlgn"
fi


#
# mapping with bowtie2 and coordinate sorting
#
# Using bowtie2
if [[ $mapper = "bowtie2" ]]
then
	echo "Using $mapper for mapping." 1>&2
	bwt2ref=`echo $ref | sed 's/\.[^\.]*$//'`
	
	bowtie2 -p $threads \
	-x $bwt2ref \
	-1 $read1 -2 $read2 \
	--rg-id $id \
	--rg SM:$sm --rg PL:illumina --rg LB:$lb --rg PU:$pu \
	$extraopts \
	| \
	samtools sort -T $outdir/samtools_tempfiles \
	-o $outdir/temp1_$outprefix.$mapper.sort.bam -;
fi

# Using bwa
if [[ $mapper = "bwa" ]]
then
	echo "Using $mapper for mapping." 1>&2
	
	bwa mem -M -t $threads \
	-R "@RG\tID:$id\tSM:$sm\tPL:illumina\tLB:$lb\tPU:$pu" \
	$extraopts \
	$ref $read1 $read2 \
	| \
	samtools sort -T $outdir/samtools_tempfiles \
	-o $outdir/temp1_$outprefix.$mapper.sort.bam -;
fi


#
# reorder reads to have the same order as the reference
#
echo "Reordering reads to match reference genome" 1>&2
java -jar $picard ReorderSam \
INPUT=$outdir/temp1_$outprefix.$mapper.sort.bam \
OUTPUT=$outdir/temp2_$outprefix.$mapper.sort.reorder.bam \
REFERENCE=$ref \
QUIET=true \
VERBOSITY=ERROR \
CREATE_INDEX=true;

if [[ $keep = "no" ]]
then
	rm $outdir/temp1*.ba*;
fi


#
# Mark duplicates with Picard and create Indel targets
#
if [[ $dups = "yes" ]]
then
	echo "Marking duplicates"  1>&2
	java -jar $picard MarkDuplicates \
	INPUT=$outdir/temp2_$outprefix.$mapper.sort.reorder.bam \
	OUTPUT=$outdir/temp3_$outprefix.$mapper.sort.reorder.markdup.bam \
	METRICS_FILE=$outdir/stats/$outprefix.bwt2.reorder.markdup_metrics \
	MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=1000 \
	CREATE_INDEX=true
	
elif [[ $dups = "no" ]]
then
	echo "Skipping marking of duplicates"  1>&2
fi

if [[ $keep = "no" ]]
then
	rm $outdir/temp2*.ba*;
fi


#
# Realign around indels
#
echo "Realigning around indels"  1>&2

if [[ $dups = "yes" ]]
then
	echo "Creating targets around Indels" 1>&2
	java -jar $gatk -T RealignerTargetCreator \
	-nt $threads -R $ref \
	-o $outdir/${output}.intervals \
	-I $$outdir/temp3_$outprefix.$mapper.sort.reorder.markdup.bam

	# Realign
	java -jar $gatk -T IndelRealigner \
	-R $ref \
	--baq CALCULATE_AS_NECESSARY \
	-I $outdir/temp3_$outprefix.$mapper.sort.reorder.markdup.bam \
	-o $outdir/${output}.bam \
	-targetIntervals $outdir/${output}.intervals \
	-noTags

elif [[ $dups = "no" ]]
then
	echo "Creating targets around Indels" 1>&2
	java -jar $gatk -T RealignerTargetCreator \
	-nt $threads -R $ref \
	-o $outdir/${output}.intervals \
	-I $$outdir/temp2_$outprefix.$mapper.sort.reorder.bam
	
	# Realign
	java -jar $gatk -T IndelRealigner \
	-R $ref \
	--baq CALCULATE_AS_NECESSARY \
	-I $outdir/temp2_$outprefix.$mapper.sort.reorder.bam \
	-o $outdir/${output}.bam \
	-targetIntervals $outdir/${output}.intervals \
	-noTags
fi

if [[ $keep = "no" ]]
then
	rm $outdir/temp3*.ba*
fi


##### Collect statistics #####
#
# collect insert size distribution
#
echo "Getting insert size distribution" 1>&2

java -jar $picard CollectInsertSizeMetrics \
INPUT=$outdir/${output}.bam \
OUTPUT=$outdir/stats/${output}.insert_size.hist \
HISTOGRAM_FILE=$outdir/stats/${output}.insert_size.pdf \
REFERENCE_SEQUENCE=$ref &


#
# samtools basic stats
#
echo "General mapping statistics (flagstat)" 1>&2

samtools flagstat $outdir/${output}.bam > $outdir/stats/${output}.flagstat &


#
# gatk callable sites stats
#
echo "Getting sites with depth > 10" 1>&2

java -jar $gatk -T CallableLoci \
-I $outdir/${output}.bam \
-R $ref \
--summary $outdir/stats/${output}.callable.summary \
-o $outdir/stats/${output}.callable.bed \
--format BED \
--minDepth 10 &

wait



