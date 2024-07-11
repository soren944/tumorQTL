
# Load necessary modules
- htslib
- singularity
- plink
- plink2
- R version 4.3.0
- python3

# Step 1: Process Gene Expression Data
- Process gene and TPM files in R
- Data process will also use the estimate method to estimate tumor purity and generate a histogram
- Requires that you already know which samples you want to use
- You can rerun if you decide to remove samples with low tumor purity
- Create a new sample list
- It will properly format gene expression data for Step 3

```
cd $home
Rscript DataProcess.R --rna_dir $rna_seq> --out_dir $home --prefix $prefix> --sample_list $samples
```

# Step 2: Compress and decompress files
```
bgzip -f JPtumor_gene_reads.gct
bgzip -f JPtumor_tpm_reads.gct
```
# Step 3: Normalize Gene Expression Data
- Uses Broad Institutes GitHub docker
- https://github.com/broadinstitute/gtex-pipeline/tree/master/qtl
- docker pull broadinstitute/gtex_eqtl:V8
  
```
singularity exec --bind $home,$home $dockerImage /bin/bash -c "/src/eqtl_prepare_expression.py ${prefix}_tpm_reads.gct.gz \
                                                              ${prefix}_gene_reads.gct.gz $GTF $samples $vcf_file ${prefix}_normalized \
                                                              --tpm_threshold 0.1 --count_threshold 6 --sample_frac_threshold 0.2 \
                                                              --normalization_method tmm"
```

# Step 4: Create position file for genes
- Necessary step to use EigenMT FDR correction
- R script create position file for the genes

```
Rscript Save_Gene_Positions.R --bed ${prefix}\_normalized.expression.bed.gz --out_dir $home --prefix $prefix --window 1e6
```

# Step 5: Process Genotype Data
- Assumes genotype data has already undergone QC
- Additional filters/thresholds
- Create genotype and genotype position files for EigenMT

```
export plink_samples=plink_sample_list.txt

for chnum in {1..22};
do
plink --bfile ${genotype_data}${chnum} --keep-allele-order --make-bed --out tmp 
plink --bfile tmp --keep $plink_samples --make-bed --maf 0.1 --biallelic-only strict --out tmp2
plink --bfile tmp2 --geno 0.02  --make-bed --out  tmp3 --biallelic-only strict
plink --bfile tmp3 --mind 0.02 --make-bed  --out tmp4 --maf 0.1
plink2 --bfile tmp4  --set-all-var-ids chr@_# --make-bed --out tmp5 --mind 0.15 --hwe 10e-8  
plink --bfile tmp5  --extract ${prefix}\_cis_position_file.txt --range --make-bed --out $prefix.chr$chnum
plink --bfile $prefix.chr$chnum --keep-allele-order --recode vcf-iid --output-chr chrMT --out $prefix.chr$chnum
Rscript Save_Genotype_Positions.R --vcf $home$prefix.chr$chnum.vcf --out_dir $home --prefix $prefix --chr $chnum 
done
```

# Step 6: Calculate PCs
- merge all genotype data to calculate genotype pcs
- identify number of pcs

export num_pcs=20

- create file which will merge genotype data from all chromosomes
```
find $home -type f -name "$prefix.chr*fam" > merge_chrom.txt
sed -i 's/\.fam$//' merge_chrom.txt
cat merge_chrom.txt
```

## merge and filter for independence
```
plink --bfile $prefix.chr1 --keep-allele-order --merge-list merge_chrom.txt --make-bed --out $prefix\_merged_full_genome --freq
plink --bfile $prefix\_merged_full_genome --indep-pairwise 50 10 0.8 --out prep_pcs_indep
plink --bfile $prefix\_merged_full_genome --extract prep_pcs_indep.prune.in --make-bed --out prep_pcs_indep2
plink --bfile prep_pcs_indep2 --pca $num_pcs --out $prefix.pcs
```
# Step 7: Calculate PEER factors
- identify number of PEER factors
- Recommendation PEER <= 10 for sample size <= 100

```
export num_peer=10
singularity exec --bind $home,$home $dockerImage2 Rscript /src/run_PEER.R $prefix\_normalized.expression.bed.gz --covariates $prefix\_tumor_impurity.txt $prefix\_PEER $num_peer -o $home
```

# Step 8: Merge and Visualize Covariates
- Covariates must be in proper format (ID column with sample IDs and the rest of the columns are the covariates)

```
export num_pcs=20
```

- Based on proportion of variance plot plot you can determine how many PCs you need
- Look for elbow
- Japanese cohort principle components did not contribute much variance

```
Rscript visualize_pca.R --out_dir $home --prefix $prefix --eigenval $prefix.pcs.eigenval --eigenvec $prefix.pcs.eigenvec --num_pca $num_pcs
```

Merge all covariates
- Can look at the generated correlation plot to determine if certain covariates have large covariance

```
Rscript merge_covariates.R --out_dir $home --prefix $prefix --pca $prefix\_pca.txt --peer $prefix\_PEER.PEER_covariates.txt --num_pca 1 --cov $cov  
```

- You must run this script again if you want to subset the covariates
- create a txt file list with correction column names including PC, PEER factor covariates
  
```
Rscript merge_covariates.R --out_dir $home --prefix $prefix --pca $prefix\_pca.txt --peer $prefix\_PEER.PEER_covariates.txt --num_pca 1 --cov $cov  --subset Y  --cov_list Covariate_names_list.txt
```
# Step 8: eQTL analysis
- This part will require submitting a job since my current code runs for long time < 24hrs for this dataset
- can specify if you want an interaction term or not and save the summary statistics for the interaction term beta

```
Rscript eQTL_analysis_pipeline_version.R -g $prefix -b $prefix\_normalized.expression.bed.gz \
                                        -p $prefix -o $home -c $prefix\_covariates_final.txt -i P
```

```
sbatch submit_eqtl_analysis.txt
```
# Step 9: FDR correction
- Need to temporarily merge all eqtl chr files to adjust p-values
- Default is Benjamini-Hochberg 
```
Rscript fdr_correction_eqtls.R --method BH  --prefix $prefix --eqtl_dir $home --out_dir $home --beta interaction 

Rscript fdr_correction_eqtls.R --method BH  --prefix $prefix --eqtl_dir $home --out_dir $home --beta genotype 
```
# Step 10: EigenMT correction
- Can download the script from https://github.com/joed3/eigenMT
- By now all our files should be formatted correctly
- I am running EigenMT Twice for each chromosome to test both the main and interaction effect
- double check you have the right python packages installed

```
for chnum in {1..22};
do
python eigenMT.py --CHROM $chnum --QTL $prefix.cis_eqtls.main_effect.chr$chnum.FDR.txt \
                                    --PHEPOS $prefix\_phe_position_file.txt \
                                    --GEN $prefix\_genotype_position_file_chr$chnum.txt \  
                                    --POS $prefix\_genotype_position_file_chr$chnum.txt \
                                    --OUT $prefix.cis_eqtls.main_effect.Eigen.chr$chnum.txt \
                                    --window 100

python eigenMT.py --CHROM $chnum --QTL $prefix.cis_eqtls.interaction_effect.chr$chnum.FDR.txt \
                                  --PHEPOS $prefix\_phe_position_file.txt \
                                  --GEN $prefix\_genotype_position_file_chr$chnum.txt \
                                  --POS $prefix\_genotype_position_file_chr$chnum.txt \
                                  --OUT $prefix.cis_eqtls.interaction_effect.Eigen.chr$chnum.txt \
                                  --window 100
done
```

# Step 11: Filter significant eQTLs, merge together, annotate
```
Rscript merge_sigpairs.R --threshold 0.05 \
                        --prefix $prefix \
                        --eqtl_dir $home --out_dir $home \
                        --beta interaction \
                        --freq ${home}${prefix}\_merged_full_genome.frq \
                        --GTF $GTF
Rscript merge_sigpairs.R --threshold 0.05 \
                        --prefix $prefix \
                        --eqtl_dir $home\
                        --out_dir $home \
                        --beta genotype \
                        --freq $home$prefix\_merged_full_genome.frq \
                        --GTF $GTF
```
# Step 12: eQTL analysis Visualization
- In Progress
```
Rscript eqtl_visualization.R --geno.beta ${home}${prefix}.sig_eqtls.main_effect.txt \
                            --int.beta ${home}${prefix}.sig_eqtls.interaction_effect.txt \
                            --out_dir $home --prefix $prefix
```
