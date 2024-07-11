#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(ggplot2)
  library(data.table)
  library(tidyverse)
  library(optparse)
})

option_list = list(
  make_option(c("-o", "--out_dir"), default=NULL, help="Directory to write files (must include backsplash at end)"),
  make_option(c("-p", "--prefix"), default=NULL, help="Prefix for all files"),
  make_option(c("-e", "--eigenvec"), default=NULL, help="plink eigenvec files"),
  make_option(c("-v", "--eigenval"), default=NULL, help="plink eigenval files"),
  make_option(c("-n", "--num_pca"), type ="numeric", default=NULL, help="number of pcs")
)

opt_parser = OptionParser(option_list=option_list)
parsed_args = parse_args(opt_parser)

# Check if eigenvec and eigenval files exist
if (!file.exists(parsed_args$eigenvec)) {
  stop(paste("Error: eigenvec file", parsed_args$eigenvec, "not found."))
}

if (!file.exists(parsed_args$eigenval)) {
  stop(paste("Error: eigenval file", parsed_args$eigenval, "not found."))
}

if (is.null(parsed_args$prefix)) {
  stop(paste("Error: prefix not found."))
}

if (is.null(parsed_args$num_pca)) {
  stop(paste("Error: num_pca not found."))
}

# Read eigenvec and eigenval files
eigenvec <- fread(parsed_args$eigenvec, header = FALSE)
eigenval <- fread(parsed_args$eigenval, header = FALSE)

# Prepare data for plotting
tmp <- eigenvec[, 2:(2 + parsed_args$num_pca)]
colnames(tmp) <- c("ID", paste0("PCA", 1:parsed_args$num_pca))

# Write PCA data to file
write.table(tmp, file = paste0(parsed_args$out_dir, parsed_args$prefix, "_pca.txt"),
            col.names = TRUE, row.names = FALSE, quote = FALSE, sep = '\t')

# Plotting
plot <- ggplot(eigenval, 
               aes(y = V1 / sum(eigenval) * 100, x = factor(1:parsed_args$num_pca), 
                   fill = factor(1:parsed_args$num_pca))) + 
  geom_bar(stat = "identity", color = "black") + 
  scale_fill_manual(values = rep(c("dodgerblue", "lightgrey"), length.out = parsed_args$num_pca)) + 
  theme_classic() + 
  labs(y = "Proportion of Variance Explained", x = "Principal Components") + 
  theme(text = element_text(size = 14, face = "bold", color = "black"),
        axis.ticks = element_line(color = "black", linewidth = 2))

# Save plot as PNG
png(filename = paste0(parsed_args$out_dir, parsed_args$prefix, "_PVE_plot.png"),
    width = 800, height = 600, res = 120, units = "px")  # Adjust dimensions and resolution as needed
print(plot)
dev.off()

# Confirm completion
print("Script completed.")