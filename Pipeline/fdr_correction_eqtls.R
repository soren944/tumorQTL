#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(optparse)
  library(tidyverse)
  library(foreach)
})

print_confirmation <- function(message) {
  cat("✅", message, "\n")
}

#note this script relies on the assumption that the SNP name has the chromosome location (i.e "chr1")
option_list = list(
  make_option(c("-d", "--eqtl_dir"), default=NULL, help="eQTL directory (must include backslash at end)"),
  make_option(c("-o", "--out_dir"), default=NULL, help="Directory to write files (must include backslash at end)"),
  make_option(c("-p", "--prefix"), default=NULL, help="Prefix for all files"),
  make_option(c("-m", "--method"), default="BH", help="R p-value correction method"),
  make_option(c("-b", "--beta"), default="genotype", help="specify (genotype) beta statistics or (interaction) statistics")
); 

opt_parser = OptionParser(option_list=option_list);
parsed_args = parse_args(opt_parser);

if (file.exists(parsed_args$eqtl_dir) && file.info(parsed_args$eqtl_dir)$isdir) {
  print_confirmation("eQTL Directory exists.")
} else {
  stop(sprintf("ERROR: eQTL Directory '%s' does not exist.", parsed_args$eqtl_dir))
}

read.eqtls <- function(x) {
  tmp <- fread(x, header = T)
  return(tmp)
}
  
if (parsed_args$beta == "genotype") {
  eqtl.files = list.files(pattern = paste0(parse_args$prefix,".cis_eqtls.main_effect.chr"), 
                          path = parsed_args$eqtl_dir)
  filename = paste0(parsed_args$out_dir,parsed_args$prefix,".cis_eqtls.main_effect.chr")
  # Function to read expression data 
} else {
  eqtl.files = list.files(pattern = paste0(parse_args$prefix,".cis_eqtls.interaction_effect.chr"), 
                                           path = parsed_args$eqtl_dir)
  filename = paste0(parsed_args$out_dir,parsed_args$prefix,".cis_eqtls.interaction_effect.chr")
}

merged.results <- map_df(eqtl.files,read.eqtls)

merged.results$FDR = p.adjust(merged.results$`p-value`,method = parsed_args$method)

for(chr_num in 1:22) {
  # Subset the data using data.table syntax and fwrite
  fwrite(merged.results[grep(paste0("chr",chr_num), merged.results$SNP)],
         file = paste0(filename, chr_num, ".FDR.txt"),
         quote = FALSE, row.names = FALSE, col.names = TRUE, sep = '\t')
}

cat("Script execution completed.")