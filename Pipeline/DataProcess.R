#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(tidyverse)
  library(data.table)
  library(estimate)
  library(ggplot2)
  library(optparse)
  library(stringr)
  library(ggplot2)
  })

print_confirmation <- function(message) {
  cat("✅", message, "\n")
}

option_list = list(
  make_option(c("-d", "--rna_dir"), default=NULL, help="RNAseq directory"),
  make_option(c("-o", "--out_dir"), default=NULL, help="Directory to write files (must include backsplash at end)"),
  make_option(c("-p", "--prefix"), default=NULL, help="Prefix for all files"),
  make_option(c("-s", "--sample_list"), default=NULL, help="List of RNAseq sample names"),
  make_option(c("-f", "--platform"), default="illumina", help="platform of RNAseq data")
); 

opt_parser = OptionParser(option_list=option_list);
parsed_args = parse_args(opt_parser);

# Confirm RNAseq directory
if (file.exists(parsed_args$rna_dir) && file.info(parsed_args$rna_dir)$isdir) {
  print_confirmation("RNAseq Directory exists.")
} else {
  stop(sprintf("ERROR: Directory '%s' does not exist.", parsed_args$rna_dir))
}

# Set RNAseq directory 
setwd(parsed_args$rna_dir)

# Function to check if file exists
check_file_exists <- function(file_path) {
  if (!file.exists(file_path)) {
    stop(sprintf("ERROR: File '%s' does not exist.", file_path))
  }
}

# Check if input files exist
check_file_exists(parsed_args$sample_list)
print_confirmation("Sample list file exists.")

# Read sample list 
samples <- fread(parsed_args$sample_list, header = TRUE)

# Function to read expression data
ReadExp <- function(x) {
  tmp <- fread(x, skip = 1)
  tmp$ID <- sub('.Aligned.sortedByCoord.out.md.bam.*',"",x)
  return(tmp)
}

# Read gene and TPM files
gene_files <- list.files(pattern = 'Aligned.sortedByCoord.out.md.bam.gene_reads.gct')
tpm_files <- list.files(pattern = 'Aligned.sortedByCoord.out.md.bam.gene_tpm.gct')

# Check IDs
missing_ids_gene <- setdiff(samples$SAMID, str_remove(gene_files, "\\.Aligned\\.sortedByCoord\\.out\\.md\\.bam\\.gene_reads\\.gct"))
missing_ids_tpm <- setdiff(samples$SAMID, str_remove(tpm_files, "\\.Aligned\\.sortedByCoord\\.out\\.md\\.bam\\.gene_tpm\\.gct"))

if (length(missing_ids_gene) > 0) {
  warning("Sample IDs not found in gene reads files:", paste(missing_ids_gene, collapse = ", "))
  stop()
}

if (length(missing_ids_tpm) > 0) {
  warning("Sample IDs not found in TPM files:", paste(missing_ids_tpm, collapse = ", "))
  stop()
}

gene_read <- map_df(gene_files, ReadExp)
tpm_read <- map_df(tpm_files, ReadExp)

# Pivot and subset data
# Write Gene Expression Matrix files
gene_read <- map_df(gene_files, ReadExp)
gene_read2 <- gene_read %>% pivot_wider(names_from = ID, values_from = Counts)
gene_final <- gene_read2[,which(colnames(gene_read2) %in% c("Name","Description",samples$SAMID))]

tpm_read <- map_df(tpm_files, ReadExp)
tpm_read2 <- tpm_read %>% pivot_wider(names_from = ID, values_from = TPM)
tpm_final <- tpm_read2[,which(colnames(tpm_read2)  %in% c("Name","Description",samples$SAMID))]

suppressWarnings({
cat("1.2",'\n',sep = '\t',file = paste0(parsed_args$out_dir,parsed_args$prefix,'_gene_reads.gct'))
cat(dim(gene_final),'\n',sep = '\t',append = T,file =  paste0(parsed_args$out_dir,parsed_args$prefix,'_gene_reads.gct'))
write.table(gene_final, paste0(parsed_args$out_dir,parsed_args$prefix,'_gene_reads.gct'),
            append = T,col.names=T,row.names=F,quote = F,sep = '\t')

cat("1.2",'\n',sep = '\t',file =  paste0(parsed_args$out_dir, parsed_args$prefix, '_tpm_reads.gct'))
cat(dim(tpm_final),'\n',sep = '\t',append = T, file =  paste0(parsed_args$out_dir, parsed_args$prefix, '_tpm_reads.gct'))
write.table(tpm_final, paste0(parsed_args$out_dir, parsed_args$prefix, '_tpm_reads.gct'),
            append = T,col.names=T,row.names=F,quote = F,sep = '\t')
})

print_confirmation("GCT files Created.")

#Remove duplicate gene symbol names and write a new gene read file for tumor purity calculation
reads.sub <- gene_final[which(!duplicated(gene_final$Description)),]
reads.sub_final <- reads.sub[,-c(1)] #remove Name column
write.table(reads.sub_final,
            file = paste0(parsed_args$out_dir,parsed_args$prefix,'_gene_reads_noduplicates.gct'),
                          col.names=T,row.names=F,quote = F,sep = '\t')
            
#Write a new gene read file filtered with common genes for tumor purity calculation  
filterCommonGenes(input.f=paste0(parsed_args$out_dir,parsed_args$prefix,'_gene_reads_noduplicates.gct'), 
                  output.f= paste0(parsed_args$out_dir, parsed_args$prefix, '_gene_reads_filtered.gct'), 
                  id="GeneSymbol")

estimateScore(input.ds = paste0(parsed_args$out_dir, parsed_args$prefix, '_gene_reads_filtered.gct'), 
              output.ds = paste0(parsed_args$out_dir, parsed_args$prefix, '_estimate_score.gct'), 
              platform = parsed_args$platform)                 

#Calculate Tumor Purity and Write out scores
scores <- fread(paste0(parsed_args$out_dir, parsed_args$prefix, '_estimate_score.gct'),header =T,skip = 2)
tumor_purity <- c("TumorPurity","TumorPurity",cos(0.6049872018+0.0001467884*as.numeric(scores[3,-c(1:2)])))
scores2 <- data.frame(SAMID = colnames(scores)[-c(1:2)], TP = tumor_purity[-c(1:2)])
scores2$SAMID <- gsub('\\.','\\-',scores2$SAMID)
print_confirmation("Estimated Tumor Purity.")

scores2$P = 1 - as.numeric(scores2$TP)
purity = as.data.frame(merge(scores2,samples,by = "SAMID"))[,c("covID","P")]

plot <- ggplot(data = scores2, aes(x = as.numeric(TP) * 100)) + 
  geom_histogram(fill = "dodgerblue", color = "black", stat = "bin") +
  theme_classic() + 
  labs(y = "Number of Samples", x = "Tumor Purity (%)") + 
  theme(
    text = element_text(size = 14, face = "bold", color = "black"),
    axis.ticks = element_line(color = "black", linewidth = 2)
  )

suppressWarnings({
png(filename = paste0(parsed_args$out_dir, parsed_args$prefix, "_tumor_purity_histogram.png"),
    width = 800, height = 600, res = 120, units = "px")  # Adjust dimensions and resolution as needed
print(plot)
dev.off()
})

write.table(purity,paste0(parsed_args$out_dir, parsed_args$prefix,"_tumor_impurity.txt"),
            col.names = T,row.names = F,quote = F, sep = '\t')

print_confirmation("Tumor Purity File and Histogram Completed.")

# Confirm completion
print_confirmation("Script execution completed.")
