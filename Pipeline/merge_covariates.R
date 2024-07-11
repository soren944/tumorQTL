#!/usr/bin/env Rscript

# Load necessary packages
suppressPackageStartupMessages({
  library(tidyverse)
  library(data.table)
  library(ggplot2)
  library(optparse)
  library(corrplot)
})

# Define command-line options
option_list <- list(
  make_option(c("-o", "--out_dir"), default=NULL, help="Directory to write files (must include backslash at end)"),
  make_option(c("-f", "--prefix"), default=NULL, help="Prefix for all files"),
  make_option(c("-p", "--peer"), default=NULL, help="Peer covariate file"),
  make_option(c("-a", "--pca"), default=NULL, help="PCA covariate file"),
  make_option(c("-s", "--subset"), default="N", help="Y or N to subset covariate file"),
  make_option(c("-l", "--cov_list"), default=NULL, help="List of final covariates to subset"),
  make_option(c("-c", "--cov"), default=NULL, help="Baseline covariates"),
  make_option(c("-n", "--num_pca"), type="integer", default=NULL, help="Number of PCA covariates")
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
check_file_exists(parsed_args$peer)
check_file_exists(parsed_args$pca)
check_file_exists(parsed_args$cov)

# Read peer, pca, and cov files
peer <- fread(parsed_args$peer, header = TRUE) %>%
  t() %>%
  as.data.frame() %>%
  rownames_to_column("ID") %>%
  rename_with(~ c("ID", "P", paste0("InferredCov", 1:(ncol(.) - 2))))

pca <- fread(parsed_args$pca, header = TRUE, select = 1:(1 + parsed_args$num_pca))
cov <- fread(parsed_args$cov, header = TRUE)

# Combine all covariates
all_cov <- reduce(list(cov, pca, peer), inner_join, by = "ID")

# Calculate correlation
corr.cov <- all_cov[, -1] %>%
  mutate_if(~ !is.numeric(.), as.numeric) %>%
  cor(method = "spearman")

# Create correlation plot
png(paste0(parsed_args$out_dir, parsed_args$prefix, '_cov_correlation_plot.png'), width = 600, height = 600, res = 120, units = "px")
corrplot(corr.cov, method = "circle", sig.level = 0.05, insig = "blank", tl.col = "black")
dev.off()

# Subset covariates if requested
if (parsed_args$subset == "Y") {
  check_file_exists(parsed_args$cov_list)
  colnames.cov <- fread(parsed_args$cov_list, header = FALSE)$V1
  write.table(as.data.frame(all_cov)[, colnames.cov],
              file = paste0(parsed_args$out_dir, parsed_args$prefix, "_covariates_final.txt"),
              col.names = TRUE, row.names = FALSE, quote = FALSE, sep = '\t')
} else {
  write.table(all_cov, file = paste0(parsed_args$out_dir, parsed_args$prefix, "_covariates_final.txt"),
              col.names = TRUE, row.names = FALSE, quote = FALSE, sep = '\t')
}

# Check for NA values in covariates
if (anyNA(all_cov)) {
  warning("There are NA values present in the covariates.")
}


