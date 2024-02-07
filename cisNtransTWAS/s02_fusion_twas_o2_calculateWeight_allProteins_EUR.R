#######################################################
## to create my own WEIGHTS for each protein-coding gene (6907 in total)

## using the pipeline built by Oliver Pain
# https://github.com/opain

## the github repo:
# https://github.com/opain/Calculating-FUSION-TWAS-weights-pipeline


rm(list = ls())
options(stringsAsFactors = FALSE)
library(reshape2)
library(dplyr)
library(data.table)
library(cowplot)
library(Biobase)
library(gdata)
## use the plasma3170-proteomics data as a test

#################################################
## 1) load info for "coordinate_file.txt"
## load in sentinel pQTLs from EUR
EUR.pQTL.plasma.block.dt <- fread('combined_2848pQTLs_addLDblockID_EUR.csv')
EUR.pQTL.plasma.block.dt
EUR.pQTL.plasma.block.dt.tidy <- EUR.pQTL.plasma.block.dt[, list(snp, Analytes, Target, chr.x, pos, type)]

### write the extraction for all features with cis and trans regions
feature.uniq.list <- unique(EUR.pQTL.plasma.block.dt.tidy$Analytes)
head(feature.uniq.list)
length(feature.uniq.list)

####################################################
## loop through all proteins for calculating the weights


# for(i in 1:10){
for(i in 1:length(feature.uniq.list)){
  cat(i, 'starts \n')
  system(paste0('/usr/bin/Rscript /home/yangc/mendelianRandomization/TWAS/fusion/washu_files/test03_plasmaProtEURcisNtrans/s1_compute_weights/d01scripts/CY_TWAS_weights_using_fusion_cisNtrans.R ',
                '--PLINK_prefix /home/yangc/associationTest/plink2_hg38/run02_plasma_proteomics_Aug2022/test10_singleFeature_mod04_callrate/plasma2338EURproteomics ',
                '--phenotype_file /home/yangc/mendelianRandomization/TWAS/fusion/washu_files/test01_plasmaProtEUR/s1_compute_weights/Calculating-FUSION-TWAS-weights-pipeline/plasmaEUR_phenotype_file.txt ',
                '--gene_name ', feature.uniq.list[i], ' ',
                '--plink /usr/local/genome/bin/plink1.9 ',
                '--gcta /home/yangc/mendelianRandomization/TWAS/fusion/fusion_twas-master/gcta_nr_robust ',
                '--gemma /home/yangc/mendelianRandomization/TWAS/fusion/gemma-0.98.5-linux-static-AMD64 ',
                '--ld_ref_dir /home/yangc/mendelianRandomization/TWAS/fusion/LDREFwashu/EUR2338 ',
                '--fusion_software /home/yangc/mendelianRandomization/TWAS/fusion/fusion_twas-master ',
                '--output_dir /home/yangc/mendelianRandomization/TWAS/fusion/washu_files/test03_plasmaProtEURcisNtrans/s1_compute_weights/d02output ',
                '--extractSNP_dir /home/yangc/mendelianRandomization/TWAS/fusion/washu_files/test03_plasmaProtEURcisNtrans/s1_compute_weights/d02output/extractSNPs'))
  
}















