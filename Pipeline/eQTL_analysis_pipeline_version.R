#!/usr/bin/env Rscript

library(data.table)
library(dplyr)
library(doParallel)
library(foreach)
library(optparse)

# Parse command line arguments
option_list <- list(
  make_option(c("-g", "--geno_prefix"), default = NULL, help = "Prefix of plink vcf files (file name format: prefix.chr1.vcf)"), 
  make_option(c("-b", "--bed"), default = NULL, help = "Normalized gene expression bed file"), 
  make_option(c("-p", "--prefix"), default = NULL, help = "Prefix"),
  make_option(c("-o", "--out_dir"), default = NULL, help = "Home and output directory must have backslash at end"),
  make_option(c("-i", "--interaction"), default = NULL, help = "Name of covariate with an interaction term with genotype"),
  make_option(c("-c", "--covariates"), default = NULL, help = "Covariate file"),
  make_option(c("-w", "--cis_window"), type = "numeric", default = 1e6, help = "length of cis window in basepairs")
)

opt_parser <- OptionParser(option_list = option_list)
parsed_args <- parse_args(opt_parser)

interaction = parsed_args$interaction
# Check if required arguments are provided
required_args <- c("geno_prefix", "bed", "covariates", "out_dir")
missing_args <- required_args[is.null(parsed_args[required_args])]
if (length(missing_args) > 0) {
  stop(sprintf("ERROR: Required argument(s) (%s) is/are missing. Please provide these arguments.", paste(missing_args, collapse = ", ")))
}

# Function to check if file exists
check_file_exists <- function(file_path) {
  if (!file.exists(file_path)) {
    stop(sprintf("ERROR: File '%s' does not exist.", file_path))
  }
}

# Check if input files exist
check_file_exists(parsed_args$bed)
check_file_exists(parsed_args$covariates)

# Read data files using fread for efficiency
gene_exp <- fread(parsed_args$bed, header = TRUE)
covar <- fread(parsed_args$covariates, header = TRUE)

# Function to perform eQTL mapping
eQTL_map <- function(chr, setgene, snp, model_data, interaction) {
  sub_model <- model_data[, c(snp, colnames(covar), "GeneExp"), with = FALSE]
  setnames(sub_model, old = snp, new = "G")
  
  if (is.null(interaction)) {
    formulaic <- formula(paste("GeneExp ~ . - ID"))
  } else {
    formulaic <- formula(paste("GeneExp ~ . + G * ", parsed_args$interaction, " - ID"))
  }
  
  model <- summary(lm(formulaic, data = sub_model))
  
  if (is.null(interaction)) {
    eqtl_results <- c(SNP = snp,
                      gene = setgene,
                      beta = signif(model$coefficients["G","Estimate"], digits = 5),
                      `t-stat` = signif(model$coefficients["G","t value"], digits = 5),
                      `p-value` = signif(model$coefficients["G","Pr(>|t|)"], digits = 5),
                      FDR = 1)
    fwrite(eqtl_results, 
           file = paste0(parsed_args$out_dir, parsed_args$prefix, ".cis_eqtls.", chr, ".txt"), 
           append = TRUE, sep = "\t")
    
  } else {
    eqtl_results.g <- data.table(SNP = snp,
                        gene = setgene,
                        beta = signif(model$coefficients["G","Estimate"], digits = 5),
                        `t-stat` = signif(model$coefficients["G","t value"], digits = 5),
                        `p-value` = signif(model$coefficients["G","Pr(>|t|)"], digits = 5),
                        FDR = 1)
    fwrite(eqtl_results.g, 
           file = paste0(parsed_args$out_dir, parsed_args$prefix, ".cis_eqtls.main_effect.", chr, ".txt"), 
           append = TRUE, sep = "\t")
    
    eqtl_results.i <- data.table(SNP = snp,
                        gene = setgene,
                        beta.int = signif(model$coefficients[paste0("G:", interaction), "Estimate"], digits = 5),
                        `t-stat.int` = signif(model$coefficients[paste0("G:", interaction), "t value"], digits = 5),
                        `p-value.int` = signif(model$coefficients[paste0("G:", interaction), "Pr(>|t|)"], digits = 5),
                        FDR = 1)
    fwrite(eqtl_results.i, 
           file = paste0(parsed_args$out_dir, parsed_args$prefix, ".cis_eqtls.interaction_effect.", chr, ".txt"), 
           append = TRUE, sep = "\t")
    
  }
}

# Parallel loop to process each chromosome
registerDoParallel(cores = detectCores())

foreach(c = 1:22, .combine = 'c', .errorhandling = "pass", .packages = c("data.table", "dplyr")) %dopar% {
  chr_num <- paste0("chr", c)
  eqtl.names = c("SNP","gene","beta","t-stat","p-value","FDR")
  if (is.null(interaction)) { 
    cat(paste(eqtl.names, sep = "\t"), '\n',
        file = paste0(parsed_args$out_dir, parsed_args$prefix, ".cis_eqtls.", chr_num, ".txt"), 
        sep = '\t', append = TRUE)
  } else {
    cat(paste(eqtl.names, sep = "\t"), '\n',
        file = paste0(parsed_args$out_dir, parsed_args$prefix, ".cis_eqtls.main_effect.", chr_num, ".txt"), 
        sep = '\t', append = TRUE)
    cat(paste(eqtl.names, sep = "\t"), '\n',
        file = paste0(parsed_args$out_dir, parsed_args$prefix, ".cis_eqtls.interaction_effect.", chr_num, ".txt"), 
        sep = '\t', append = TRUE)
  }
  # Read genotype data for the current chromosome using fread
  genotype_chr <- fread(paste0(parsed_args$out_dir, parsed_args$geno_prefix, ".", chr_num, ".vcf"), 
                        header = TRUE, drop = c("QUAL", "FILTER", "INFO", "FORMAT", "REF", "ALT"))
  
  # Replace any colons and dashes from variant name since it causes issues later
  if (any(grepl("[-:]", genotype_chr$ID))) {
    genotype_chr$ID <- gsub("[-:]", "_", genotype_chr$ID)
  }
  
  # Subset gene expression data for the current chromosome
  Gene_exp_chr <- gene_exp[`#chr` == chr_num, ]
  
  # Convert genotype positions to numeric and set row names
  genotype_chr$POS <- as.numeric(genotype_chr$POS)
  row.names(genotype_chr) <- genotype_chr$ID
  
  # Log total number of genes for the current chromosome
  cat(paste("Total # of Genes:", nrow(Gene_exp_chr), sep = " "), '\n', 
      file = paste0(parsed_args$out_dir, parsed_args$prefix, ".cis_eqtls.", chr_num, ".log"))
  
  # Iterate over genes in the current chromosome
  foreach(gene = 1:nrow(Gene_exp_chr), .combine = 'c', .errorhandling = "pass") %do% {
    setgene <- Gene_exp_chr$gene_id[gene]
    start <- Gene_exp_chr$start[gene]  
    end <- Gene_exp_chr$end[gene]
    
    # Subset genotype data around the gene's position
    CisGeno <- genotype_chr[POS >= start - parsed_args$cis_window & POS <= end + parsed_args$cis_window, ]
    snplist <- CisGeno$ID
    
    if (nrow(CisGeno) >= 1) {
      # Transpose the subsetted data.table and assign SNP names to ID column
      CisGeno_transposed <- transpose(CisGeno[, -c(1:2)])
      colnames(CisGeno_transposed) <- unlist(CisGeno_transposed[1, ])
      CisGeno_transposed <- CisGeno_transposed[-1, ]
      CisGeno_transposed[, ID := colnames(CisGeno)[-c(1:3)]]
      
      # Convert genotype columns
      convert_genotype <- function(x) {
        x <- as.character(x)
        x[x == "0/0"] <- 0
        x[x == "0/1"] <- 1
        x[x == "1/1"] <- 2
        x[x == "./."] <- NA
        as.numeric(x)   }
      
      genotype_cols <- setdiff(names(CisGeno_transposed), c("ID"))
      # Apply convert_genotype function to genotype columns
      CisGeno_transposed[, (genotype_cols) := lapply(.SD, function(x) {
        x <- convert_genotype(x)
        mean_val <- mean(x, na.rm = TRUE)
        x[is.na(x)] <- mean_val
        x
      }), .SDcols = genotype_cols]
      
      # Create a combined dataset for modeling
      GeneExpSub <- data.table(ID = colnames(gene_exp)[5:ncol(gene_exp)], 
                               GeneExp = unlist(Gene_exp_chr[gene, .SD, .SDcols = -c(1:4)]))
      model_data <- Reduce(function(x, y) merge(x, y, by = "ID", all = TRUE), 
                           list(GeneExpSub, CisGeno_transposed, covar))
      
      # Create a combined dataset for modeling
      GeneExpSub <- data.table(ID = colnames(gene_exp)[5:ncol(gene_exp)], 
                               GeneExp = unlist(Gene_exp_chr[gene, .SD, .SDcols = -c(1:4)]))
      model_data <- Reduce(function(x, y) merge(x, y, by = "ID", all = TRUE), 
                           list(GeneExpSub, CisGeno_transposed, covar))
      
      # Perform eQTL mapping for each SNP
      foreach(snp = snplist, .combine = 'c', .errorhandling = "pass") %do% {
        eQTL_map(chr_num, setgene, snp, model_data, parsed_args$interaction)
      }
    }
    # Log completion of processing for the current gene
    cat(paste("Finished Gene:", gene, sep = " "), '\n', append = TRUE, 
        file = paste0(parsed_args$out_dir, parsed_args$prefix, ".cis_eqtls.", chr_num, ".log"))
  }
}
