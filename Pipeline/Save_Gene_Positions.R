#!/usr/bin/env Rscript
.libPaths('~/R/x86_64-pc-linux-gnu-library/4.3/')
library(data.table)
library(stringr)
library(optparse)

option_list = list(
  make_option(c("-b", "--bed"), default=NULL, help="normalized gene expression bed file"),
  make_option(c("-o", "--out_dir"), default=NULL, help="Directory to write files (must include backsplash at end)"),
  make_option(c("-p", "--prefix"), default=NULL, help="Prefix for all files"),
  make_option(c("-w", "--window"), default=1e6, help="cis-window length")
); 

opt_parser = OptionParser(option_list=option_list);
parsed_args = parse_args(opt_parser);
cis = as.numeric(parsed_args$window)

gx <- fread(parsed_args$bed,header = T)
gx$`#chr` = gsub("chr","", gx$`#chr`)

tmp = as.data.frame(gx)[,c("#chr","start","end","gene_id")]
tmp$end = format(tmp$end + cis, scientific = FALSE,nsmall = 0)
tmp$start = format(pmax(tmp$start - cis,0), scientific = FALSE,nsmall = 0)

write.table(tmp, paste0(parsed_args$out_dir,parsed_args$prefix,"_cis_position_file.txt"),
            col.names = F, row.names = F,sep = "\t",quote = F)

tmp2 = as.data.frame(gx)[,c("gene_id","#chr","start","end")]
colnames(tmp2) = c("gene_id","chrom_probe","s1","s2")
write.table(tmp2,paste0(parsed_args$out_dir,parsed_args$prefix,"_phe_position_file.txt"),col.names = F, row.names = F,sep = "\t",quote = F)




