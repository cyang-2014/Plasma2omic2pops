#######################################################
## to create my own WEIGHTS for each protein-coding gene (6907 in total)

## using the pipeline built by Oliver Pain
# https://github.com/opain

## the github repo:
# https://github.com/opain/Calculating-FUSION-TWAS-weights-pipeline

## update OP_TWAS_weights_using_fusion.R into CY_TWAS_weights_using_fusion_cisNtrans.R
rm(list = ls())
options(stringsAsFactors = FALSE)
library(reshape2)
library(dplyr)
library(data.table)
library(cowplot)
library(Biobase)
library(gdata)
library(foreach)
## use the plasma3170-proteomics data as a test
###########################################################################

################################################################################################
#### to extract all cis and trans region given the same feature for getting one tmp.bfiles
## load in sentinel pQTLs from EUR
EUR.pQTL.plasma.block.dt <- fread('combined_2848pQTLs_addLDblockID_EUR.csv')
EUR.pQTL.plasma.block.dt
EUR.pQTL.plasma.block.dt.tidy <- EUR.pQTL.plasma.block.dt[, list(snp, Analytes, Target, chr.x, pos, type)]




## load the full genome-wide variant lists
EUR.pQTL.variants.dt <- fread('/home/yangc/associationTest/plink2_hg38/run02_plasma_proteomics_Aug2022/test10_singleFeature_mod04_callrate/plasma2338EURproteomics.bim')
colnames(EUR.pQTL.variants.dt) <- c('chr', 'variantID', 'cM', 'pos', 'A1', 'A2')

DIRout_extractSNP <- '/home/yangc/mendelianRandomization/TWAS/fusion/washu_files/test03_plasmaProtEURcisNtrans/s1_compute_weights/d02output/extractSNPs/'

## get the variant list for cis and trans regions per feature, e.g. proteinID X16300.4
singleFeature.dt.regions <- EUR.pQTL.plasma.block.dt.tidy[Analytes == 'X16300.4']
singleFeature.dt.regions
singleFeature.dt_regions.extract <- foreach(i = 1:nrow(singleFeature.dt.regions), .combine = 'rbind') %do% {
  singleFeature.dt_regionI <- singleFeature.dt.regions[i]
  EUR.pQTL.variants.dt.sub <- EUR.pQTL.variants.dt[chr == singleFeature.dt_regionI$chr.x][pos > singleFeature.dt_regionI$pos-1e6 &
                                                                                            pos < singleFeature.dt_regionI$pos+1e6]
  return(EUR.pQTL.variants.dt.sub)
  
}

fwrite(singleFeature.dt_regions.extract[, list(variantID)],
  paste0(DIRout_extractSNP,
         'X16300.4', '_EUR.txt'), col.names = FALSE)



# Extract gene from phenotype file
system(paste('mkdir ',opt$output_dir,'/temp',sep=''))
system(paste('awk -f /home/yangc/mendelianRandomization/TWAS/fusion/washu_files/test03_plasmaProtEURcisNtrans/s1_compute_weights/d01scripts/t.awk c1=FID c2=IID c3=',
             'X16300.4',' ',
             '/home/yangc/mendelianRandomization/TWAS/fusion/washu_files/test01_plasmaProtEUR/s1_compute_weights/Calculating-FUSION-TWAS-weights-pipeline/plasmaEUR_phenotype_file.txt',' > ',
             '/home/yangc/mendelianRandomization/TWAS/fusion/washu_files/test03_plasmaProtEURcisNtrans/s1_compute_weights/d02output',
             '/temp/temp_', 'X16300.4','.pheno', sep=''))
# Using PLINK, extract variants for the cis+trans region of each feature

system(paste0('/usr/local/genome/bin/plink1.9 ',
              '--bfile /home/yangc/associationTest/plink2_hg38/run02_plasma_proteomics_Aug2022/test10_singleFeature_mod04_callrate/plasma2338EURproteomics ',
              '--make-bed --pheno ',
              '/home/yangc/mendelianRandomization/TWAS/fusion/washu_files/test03_plasmaProtEURcisNtrans/s1_compute_weights/d02output',
              '/temp/temp_', 'X16300.4','.pheno',
              ' --geno 0.02 ', '--extract ',
              DIRout_extractSNP, 'X16300.4', '_EUR.txt', sep=''))

###################################################
### write the extraction for all features with cis and trans regions
feature.uniq.list <- unique(EUR.pQTL.plasma.block.dt.tidy$Analytes)
head(feature.uniq.list)
length(feature.uniq.list)

for(feature_j in feature.uniq.list) {
   singleFeature.dt.regions <- EUR.pQTL.plasma.block.dt.tidy[Analytes == feature_j]
   singleFeature.dt.regions
   singleFeature.dt_regions.extract <- foreach(i = 1:nrow(singleFeature.dt.regions), .combine = 'rbind') %do% {
     singleFeature.dt_regionI <- singleFeature.dt.regions[i]
     EUR.pQTL.variants.dt.sub <- EUR.pQTL.variants.dt[chr == singleFeature.dt_regionI$chr.x][pos > singleFeature.dt_regionI$pos-1e6 &
                                                                                               pos < singleFeature.dt_regionI$pos+1e6]
     return(EUR.pQTL.variants.dt.sub)
     
   }
   
   fwrite(singleFeature.dt_regions.extract[, list(variantID)],
          paste0(DIRout_extractSNP,
                 feature_j, '_EUR.txt'), col.names = FALSE)
 }



