#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(data.table)
  library(optparse)
  library(tidyverse)
})

print_confirmation <- function(message) {
  cat("✅", message, "\n")
}

#note this script relies on the assumption that the SNP name has the chromosome location (i.e "chr1")
option_list = list(
  make_option(c("-d", "--eqtl_dir"), default=NULL, help="eQTL directory (must include backslash at end)"),
  make_option(c("-o", "--out_dir"), default=NULL, help="Directory to write files (must include backslash at end)"),
  make_option(c("-p", "--prefix"), default=NULL, help="Prefix for all files"),
  make_option(c("-t", "--threshold"), default=0.05,type = "numeric", help="significance threshold"),
  make_option(c("-f", "--freq"), default=NULL, help="Merged genotype plink .frq file"),
  make_option(c("-g", "--GTF"), default=NULL, help="Gencode GTF file"),
  make_option(c("-b", "--beta"), default="genotype", help="specify (genotype) beta statistics or (interaction) statistics")
); 

opt_parser = OptionParser(option_list=option_list);
parsed_args = parse_args(opt_parser);

if (file.exists(parsed_args$eqtl_dir) && file.info(parsed_args$eqtl_dir)$isdir) {
  print_confirmation("eQTL Directory exists.")
} else {
  stop(sprintf("ERROR: eQTL Directory '%s' does not exist.", parsed_args$eqtl_dir))
}

gtf <- fread(parsed_args$GTF,header = F,select = c("V3","V4","V5","V7","V9"))
colnames(gtf)[1] <- "type"
gtf <- subset(gtf,type == "gene")
colnames(gtf) <- c("type","TSS","TES","strand","info")
gtf_info <- str_split_fixed(gsub('\"',"",gtf_sub$info),"; ",10)
gtf_info <- as.data.frame(gtf_info)[,c(1,3,4)]
colnames(gtf_info) <- c("gene","gene_type","gene_name")
gtf_info$gene_id <- gsub("gene ","",gtf_info$gene)
gtf_info$gene_type <- gsub("gene_type ","",gtf_info$gene_type)
gtf_info$gene_name <- gsub("gene_name ","",gtf_info$gene_name)
gtf_full <- cbind(gtf[,!c("info")],gtf_info)

freq <- fread(parsed_args$freq,header = T,select = c(1:5))
colnames(freq) <- c("chr","SNP","alt_allele","ref_allele","maf")

read.eqtls <- function(x) {
  tmp <- fread(x, header = T)
  return(tmp)
}
  
if (parsed_args$beta == "genotype") {
    eqtl.files = list.files(pattern = paste0(parse_args$prefix,".cis_eqtls.main_effect.Eigen.chr"), 
                            path = parsed_args$eqtl_dir)
    filename = paste0(parsed_args$out_dir,parsed_args$prefix,".sig_eqtls.main_effect.txt")
    # Function to read expression data 
} else {
    eqtl.files = list.files(pattern = paste0(parse_args$prefix,".cis_eqtls.interaction_effect.Eigen.chr"), 
                            path = parsed_args$eqtl_dir)
    filename = paste0(parsed_args$out_dir,parsed_args$prefix,".sig_eqtls.interaction_effect.txt")
}

merged.results <- map_df(eqtl.files,read.eqtls)

sig.res <- merged.results %>% filter(BF <= parsed_args$threshold)

sig.res.annotated = reduce(list(sig.res, freq, gtf_full), inner_join )

fwrite(sig.res.annotated, file = filename,
         quote = FALSE, row.names = FALSE, col.names = TRUE, sep = '\t')

cat("Script execution completed.")


