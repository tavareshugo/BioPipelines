suppressMessages({
  library(optparse)
})

#
# Define options ----
#
option_list = list(
  make_option(c("-b", "--bam"), type = "character", default = NULL, 
              help = "The alignment .bam file", metavar = "file"),
  make_option(c("-o", "--out"), type="character", default = NULL, 
              help = "Output file name (the output format will be a CSV file). 
              If not given, will output to the same directory as the input file 
              with extension '.cts.csv'", 
              metavar = "file"),
  make_option("--gtf", type = "character", default = NULL,
              help = "The annotation file in gtf format. If you have a gff3 
              file convert it to gtf first. For example, you can do this using 
              the 'gffread' function from cufflinks, like so: 
                gffread annotation.gff3 -T -o annotation.gtf",
              metavar = "character"),
  make_option("--mode", type = "character", default = "Union",
              help = "The method for counting reads overlapping with the exons. 
              One of (case-sensitive): Union, IntersectionStrict, IntersectionNotEmpty. 

              See more details on Figure 1 of the GenomicAligments vignette: 
              https://bioconductor.org/packages/devel/bioc/vignettes/GenomicAlignments/inst/doc/summarizeOverlaps.pdf", 
              metavar = "character"),
  make_option("--single", action = "store_true", default = FALSE,
              help = "Use this option if your library is single-end."),
  make_option("--stranded", action = "store_true", default = FALSE,
              help = "Use this option if your library is strand-specific."),
  make_option("--fragments", action = "store_true", default = FALSE,
              help = "Use this option if you want to count same-strand pairs, 
              singletons, reads with unmapped pairs and other fragments. 
              For more details see the 'fragments' option in the summarizeOverlaps 
              function of the GenomicAlignments package.")
)

opt_parser = OptionParser(option_list=option_list,
                          description = "
This script counts the number of reads overlapping the 
exons from an annotation file. It uses the 'summariseOverlaps' 
function from the Bioconductor package GenomicAlignments.")
opt = parse_args(opt_parser)


# # For testing purposes
# opt <- list(bam = "~/mount/hpc_slcu/temp/test_rnaseqcount/test_alignment.bam",
#             gtf = "~/mount/hpc_slcu/reference/arabidopsis/tair10/annotation/Arabidopsis_thaliana.TAIR10.34.gtf",
#             mode = "Union")

#
# Check options ----
#
# Check that the mode option is valid
if(!(opt$mode %in% c("Union", "IntersectionStrict", "IntersectionNotEmpty"))){
  stop("The count method (--mode option) has to be one of Union, 
       IntersectionStrict or IntersectionNotEmpty. ", opt$mode, " is not valid.")
}

# Check that the bam input exists
if(is.null(opt$bam)){
  print_help(opt_parser)
  stop("Please provide an input bam file.")
} else if(!file.exists(opt$bam)){
  stop("Bam file could not be found: ", opt$bam)
}

# Check if the output file exists (or directory it should be placed on)
if(is.null(opt$out)){
  # Make output file same as input file with new extension
  opt$out <- paste0(tools::file_path_sans_ext(opt$bam), ".cts.csv")
  
  # Give a warning
  warning("No output file provided. Will save the file as: ", opt$out)
  
} else if(!dir.exists(dirname(opt$out))){
  stop("Output directory does not exist: ", dirname(opt$out))
}


# Check annotation file
if(is.null(opt$gtf)){
  print_help(opt_parser)
  stop("Please provide an input gtf annotation file.")
} else if(!file.exists(opt$gtf)){
  stop("GTF file could not be found: ", opt$gtf)
} else if(tools::file_ext(opt$gtf) != "gtf"){
  stop("The annotation file has to be GTF format. Your file's extension is: ", 
       tools::file_ext(opt$gtf))
}


#
# Read counts ----
#
suppressMessages({
  library(GenomicAlignments)
  library(Rsamtools)
  library(GenomicFeatures)
})


cat("Reading annotation file...\n")
# Read annotation file
annotation <- makeTxDbFromGFF(opt$gtf, format = "gtf")

# Extract exons grouped by gene
exons <- exonsBy(annotation, by="gene")

# Print some information
cat("There are", summary(exons)[1], "genes in your annotation.\n")

# Create reference to BAM file
bam_file <- BamFile(opt$bam)

# Count reads
cat("Counting reads...\n")
counts <- summarizeOverlaps(features = exons, 
                            reads = bam_file,
                            mode = opt$mode,
                            singleEnd = opt$single,
                            ignore.strand = !opt$stranded,
                            fragments = opt$fragments)


#
# Prepare output table ----
#
suppressMessages(library(dplyr))

# Make a data frame with gene features
exonsdf <- unlist(exons)
exonsdf <- data.frame(seqnames = seqnames(exonsdf),
                      start = start(exonsdf)-1,
                      end = end(exonsdf),
                      gene_id = names(exonsdf))

# Summarise by gene (with start and end positions)
genesdf <- exonsdf %>% 
  group_by(gene_id, seqnames) %>%
  summarise(gene_start = min(c(start, end)),
            gene_end = max(c(start, end)),
            gene_length = (gene_end - gene_start + 1)/1000)


# Make a data.frame of read counts
countsdf <- data.frame(gene_id = rownames(counts),
                       counts = as.integer(assay(counts)))

# Join with the gene table
countsdf <- left_join(countsdf, genesdf, by = "gene_id")

# Calculate FPKM and TPM
countsdf <- countsdf %>% 
  mutate(rpm = counts/(sum(counts)/1e6),
         fpkm = rpm/gene_length,
         rpk = counts/gene_length,
         tpm = rpk/(sum(rpk)/1e6)) %>% 
  select(-rpm, -rpk)

# Write output
write.csv(countsdf, opt$out, row.names = FALSE)


