#!/usr/bin/Rscript

suppressMessages({
  if(!require(data.table)){
    stop("please install the data.table package: install.packages(\"data.table\"")
  }
  
  if(!require(optparse)){
    stop("please install the `optparse` package: install.packages(\"optparse\")")
  }
})


#
# Define options ----
#
option_list = list(
  make_option(c("--snape"), type = "character", default = NULL, 
              help = "A file produced by SNAPE.", metavar = "file"),
  make_option("--snp_list", type = "character", default = NULL,
              help = "A file with a list of known SNPs. This file should have two columns: 
              chromossome name and position. The format is automatically detected.",
              metavar = "file"),
  make_option("--snp_number", type = "integer", default = 10000L,
              help = "An integer of the number of SNPs to subsample the original file.
              This option is ignored if --snp_list is provided. If the number is greater 
              than the total number of SNPs in the input file, it is also ignored and 
              all sites will be written.",
              metavar = "file"),
  make_option("--seed", type = "integer", default = NULL,
              help = "An integer for setting the seed for random number generator,
              in case reproducibility is important."),
  make_option(c("--out"), type="character", default = NULL, 
              help = "Output file name (the output format will be a CSV file). 
              If not given, will output to the same directory as the input file 
              with extension '.snape.csv'", 
              metavar = "file")
  )

opt_parser = OptionParser(option_list=option_list,
                          description = "
                          This script reads and filters snape-pooled output to include 
                          either a subset of user-provided SNPs or a random subset. The 
                          output is a CSV format with appropriate column names.")
opt = parse_args(opt_parser)


#
# Check options ----
#
# Check if input file exists
if(is.null(opt$snape)){
  print_help(opt_parser)
  stop("Please provide an input snape-pooled file.")
} else if(!file.exists(opt$snape)){
  print_help(opt_parser)
  stop("snape-pooled file could not be found: ", opt$snape)
}

# Check if output file directory exists
if(is.null(opt$out)){
  # Make output file same as input file with new extension
  opt$out <- paste0(tools::file_path_sans_ext(opt$snape), ".snape.csv")
  
  # Give a warning
  warning("No output file provided. Will save the file as: ", opt$out)
  
} else if(!dir.exists(dirname(opt$out))){
  stop("Output directory does not exist: ", dirname(opt$out))
}

# Check SNP subsetting options
if(is.null(opt$snp_list)){
  if(!is.integer(opt$snp_number)){
    print_help(opt_parser)
    stop("--snp_number has to be an integer number. It is: ", opt$snp_number)
  }
  
  message("Randomly selecting ", opt$snp_number, " SNPs from input file.")
  
} else if(!file.exists(opt$snp_list)){
  stop("SNP list file could not be found: ", opt$snp_list)
}


#
# Read SNP file ----
#
message("Reading snape-pooled file: ", opt$snape)

snape <- fread(opt$snape,
               col.names = c("chrom", "pos", "ref_allele", "ref_count", "alt_count", 
                             "ref_qual", "alt_qual", "alleles", "prob_polymorphic", 
                             "prob_fixed", "alt_freq"),
               colClasses = c("character", "integer", "character", "integer", "integer",
                              "integer", "integer", "character", "numeric", "numeric", "numeric"),
               na.strings = "*",
               showProgress = FALSE,
               integer64 = "numeric")



#
# Subset SNPs and write output ----
#
message("Subsetting SNPs...")

# Subset SNPs
if(!is.null(opt$snp_list)){  # If SNP list is provided
  # Read SNP list
  snp_list <- fread(opt$snp_list, showProgress = FALSE)
  
  # Get SNP ID
  snp_subset <- paste(snp_list[[1]], snp_list[[2]], sep = "_")
  
  # Subset original table
  snape <- snape[paste(chrom, pos, sep = "_") %in% snp_subset, ]
  
  if(nrow(snape) == 0){
    stop("No SNPs left after filtering. 
         Please check that your SNP list file is in the right format: 
         first column chromosome, second column position")
  } 
  
} else if(opt$snp_number >= nrow(snape)){  # If number of SNPs requested is greater than total
  
  message("The number of SNPs provided with --snp_number is greater than the total 
          number of SNPs in the input file. Writting all SNPs to output.")
  
} else {  # If number of SNPs requested is OK
  
  set.seed(opt$seed)
  snape <- snape[sample(nrow(snape), opt$snp_number), ]
  
}

# Write output
message("Writting output file: ", opt$out)
write.csv(snape, file = opt$out, row.names = FALSE)
