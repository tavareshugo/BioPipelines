# BioPipelines

Personal scripts to run general and often-used pipelines.
In all cases, it might be necessary to adjust variables with paths to 
GATK and Picard Tools.


### index_reference.sh

Indexes a reference genome with `bowtie2-build`, `samtools faidx`, 
`CreateSequenceDictionary` from Picard Tools and `bwa index`.

Instructions of usage can be obtained by running the script with no 
options: `bash index_reference.sh`


### mapping_pipeline.sh

Does mapping and post-mapping steps using several tools, including 
indel-realignment.

Instructions of usage can be obtained by:
`bash mapping_pipeline.sh -h`


### qc_fastq.sh

Quality control and filtering using FastQC and cutadapt. Uses GNU 
parallel to parallelize filtering of the reads.

Instructions of usage can be obtained by:
`bash qc_fastq.sh -h`


### count_overlaps_bam_gtf.R

An R script to count read alignments overlapping an annotation (GTF) file.

It does so by using the 
[GenomicAlignments](https://bioconductor.org/packages/release/bioc/html/GenomicAlignments.html), 
[Rsamtools](https://bioconductor.org/packages/release/bioc/html/Rsamtools.html), 
and [GenomicFeatures](https://bioconductor.org/packages/release/bioc/html/GenomicFeatures.html) 
packages. 

The pipeline is based on [this BioConductor tutorial](https://www.bioconductor.org/help/workflows/rnaseqGene/).

Instructions for usage can be obtained by:
`Rscript --help`
