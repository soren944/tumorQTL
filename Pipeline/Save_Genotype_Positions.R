#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(data.table)
  library(stringr)
  library(optparse)
})


# Define command-line options
option_list <- list(
  make_option(c("-v", "--vcf"), default=NULL, help="Plink vcf file"),
  make_option(c("-o", "--out_dir"), default=NULL, help="Directory to write files (must include backslash at end)"),
  make_option(c("-p", "--prefix"), default=NULL, help="Prefix for all files"),
  make_option(c("-c", "--chr"), default=NULL, help="Specify chromosome for genotype position file")
)

# Parse command-line arguments
opt_parser <- OptionParser(option_list=option_list)
parsed_args <- parse_args(opt_parser)

# Function to check if file exists
check_file_exists <- function(file_path) {
  if (!file.exists(file_path)) {
    stop(sprintf("ERROR: File '%s' does not exist.", file_path))
  }
}

# Check if input files exist
check_file_exists(parsed_args$vcf)

# Read genotype data
geno <- fread(parsed_args$vcf, header = TRUE, drop = c("QUAL", "FILTER", "INFO", "FORMAT", "REF", "ALT"))

# Extract SNP information for the specified chromosome
tmp <- geno[, .(snp = ID, chr_snp = `#CHROM`, pos = POS)]

# Write output file for SNP information
output_file <- paste0(parsed_args$out_dir, parsed_args$prefix, "_genotype_position_file_chr", parsed_args$chr, ".txt")
fwrite(tmp, file=output_file, col.names=TRUE, row.names=FALSE, sep="\t", quote=FALSE)

# Define function to convert genotype
convert_genotype <- function(x) {
  x <- as.character(x)
  x[x == "0/0"] <- 0
  x[x == "0/1"] <- 1
  x[x == "1/1"] <- 2
  x[x == "./."] <- NA
  as.numeric(x)
}

# Apply convert_genotype function to genotype columns
suppressWarnings({
genotype_cols <- setdiff(names(geno), "ID")
geno[, (genotype_cols) := lapply(.SD, convert_genotype), .SDcols = genotype_cols]})

# Write output file for genotype data
output_file2 <- paste0(parsed_args$out_dir, parsed_args$prefix, "_genotype_file_chr", parsed_args$chr, ".txt")
fwrite(geno, file=output_file2, col.names=TRUE, row.names=FALSE, sep="\t", quote=FALSE)

# Print confirmation message
cat(sprintf("Genotype position file for chromosome %s saved to: %s\n", parsed_args$chr, output_file))
cat(sprintf("Genotype file for chromosome %s saved to: %s\n", parsed_args$chr, output_file2))