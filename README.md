#bioPipelines
Personal shell scripts to run general and often-used pipelines.
In all cases, it might be necessary to adjust variables with paths to GATK and Picard Tools.

###index_reference.sh
Indexes a reference genome with `bowtie2-build`, `samtools faidx`, `CreateSequenceDictionary` from Picard Tools and `bwa index`.

Instructions of usage can be obtained by running the script with no options:
`bash index_reference.sh`


###mapping_pipeline.sh
Does mapping and post-mapping steps using several tools, including indel-realignment.

Instructions of usage can be obtained by:
`bash mapping_pipeline.sh -h`


