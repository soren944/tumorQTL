#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(data.table)
  library(optparse)
  library(tidyverse)
  library(ggplot2)
})


option_list = list(
  make_option(c("-e", "--geno.beta"), default=NULL, help="main effect eQTL annotated results"),
  make_option(c("-e", "--int.beta"), default=NULL, help="interaction effect eQTL annotated results"),
  make_option(c("-o", "--out_dir"), default=NULL, help="Directory to write files (must include backslash at end)"),
  make_option(c("-p", "--prefix"), default=NULL, help="Prefix for all files")
)

opt_parser = OptionParser(option_list=option_list);
parsed_args = parse_args(opt_parser);

if (file.exists(parsed_args$eqtl_dir) && file.info(parsed_args$eqtl_dir)$isdir) {
  print_confirmation("eQTL Directory exists.")
} else {
  stop(sprintf("ERROR: eQTL Directory '%s' does not exist.", parsed_args$eqtl_dir))
}

sig.int = fread(parsed_args$geno.beta,header = T)
sig.geno = fread(parsed_args$int.beta,header = T)
overlap.sig = merge(sig.geno,sig.int,
                by = c("SNP","gene","gene_name","gene_type","maf","ref_allele",
                       "alt_allele", "chr"), suffixes = c("",".int"))

sig.int$chr <- factor(sig.int$chr,levels = c(1:22))
sig.geno$chr <- factor(sig.geno$chr,levels = c(1:22))


cat("Number of cis-eQTLs/eGenes:",dim(sig.geno)[1],"\n")
cat("Number of signficant interaction effects:",dim(sig.int)[1],"\n")
cat("Overlap: eQTLs with significant interaction effects:",dim(overlap.sig)[1],"\n")

sig.geno$interaction = "p-value > 0.05"
sig.geno$interaction[which(sig.geno$gene %in% overlap.sig$gene)] = "p-value <= 0.05"


fig1 = ggplot(sig.geno,aes(x = chr,y = ..count..,fill = interaction)) + geom_bar(show.legend = T,color = "black") +
  scale_fill_manual( values = c("dodgerblue","lightgrey") ) + 
  theme_classic() + labs(x = "Chromosome", y = "Number of eQTLs",fill = "Interaction Term") +
  theme(text = element_text(size = 14, face = "bold", color = "black"),
        axis.ticks = element_line(color = "black", linewidth = 2))

suppressWarnings({
  png(filename = paste0(parsed_args$out_dir, parsed_args$prefix, "_eqtl_count_by_chr_barplot.png"),
      width = 800, height = 500, res = 120, units = "px")  # Adjust dimensions and resolution as needed
  print(fig1)
  dev.off()
})


fig2 = ggplot(sig.geno,aes(x = gene_type,y = log10(..count..),fill = interaction)) + 
  geom_bar(show.legend = T,color = "black") +
  scale_fill_manual( values = c("dodgerblue","lightgrey") ) + 
  theme_classic() + labs(x = "Gene Type", y = "log10(# cis-eQTL)",fill = "Interaction Term")  +
  theme(axis.text.x = element_text(angle = 90, vjust = 1, hjust=1)) +
  coord_flip() + 
  theme(text = element_text(size = 14, face = "bold", color = "black"),
        axis.ticks = element_line(color = "black", linewidth = 2))

suppressWarnings({
  png(filename = paste0(parsed_args$out_dir, parsed_args$prefix, "_eqtl_type_by_chr_barplot.png"),
      width = 900, height = 600, res = 120, units = "px")  # Adjust dimensions and resolution as needed
  print(fig2)
  dev.off()
})


fig3 <- ggplot(data = sig.geno, aes(x = beta)) + 
  geom_histogram(fill = "dodgerblue", color = "black",) +
  theme_classic() + 
  labs(y = "eQTL Main Effect Size", x = "Frequency") + 
  theme(
    text = element_text(size = 14, face = "bold", color = "black"),
    axis.ticks = element_line(color = "black", linewidth = 2)
  )

suppressWarnings({
  png(filename = paste0(parsed_args$out_dir, parsed_args$prefix, "_eQTL_histogram_effect_sizes.png"),
      width = 800, height = 600, res = 120, units = "px")  # Adjust dimensions and resolution as needed
  print(plot)
  dev.off()
  })
