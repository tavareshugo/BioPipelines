#bioPipelines
General shell scripts to run often-used pipelines.

###index_reference.sh
Indexes a reference genome with `bowtie2-build`, `samtools faidx`, `CreateSequenceDictionary` from Picard Tools and `bwa index`.

Instructions of usage can be obtained by running the script with no options:
`bash index_reference.sh`


###mapping_pipeline.sh
Does mapping and post-mapping steps using several tools, including indel-realignment.

Instructions of usage can be obtained by:
`bash mapping_pipeline.sh -h`


