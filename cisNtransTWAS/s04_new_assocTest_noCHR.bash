#!/bin/bash

OUTdir=/home/yangc/mendelianRandomization/TWAS/fusion/washu_files/test03_plasmaProtEURcisNtrans/s2_perform_expression_imputation/check03eurT2D2022_fiziIMPUTEsumstats
LDrefDIR=/home/yangc/associationTest/plink2_hg38/run02_plasma_proteomics_Aug2022/test10_singleFeature_mod04_callrate
weightDIR=/home/yangc/mendelianRandomization/TWAS/fusion/washu_files/test03_plasmaProtEURcisNtrans/s1_compute_weights/d02output/TWAS_weights_package
weightTEMPdir=/home/yangc/mendelianRandomization/TWAS/fusion/washu_files/test03_plasmaProtEURcisNtrans/s1_compute_weights/d02output/temp
DISEASEpath=/home/yangc/mendelianRandomization/TWAS/fizi/output/combined_new_imputed.cleaned_eurT2D2022.tsv.gz


Rscript /home/yangc/mendelianRandomization/TWAS/fusion/fusion_twas-master/CY.FUSION.assoc_test.R \
--sumstats $DISEASEpath \
--weights $weightDIR/TWAS_Weights.pos \
--weights_dir $weightDIR \
--weightTEMPdir $weightTEMPdir \
--ref_ld $LDrefDIR/plasma2338EURproteomics \
--out $OUTdir/output/NEW4_eurT2Drisk2022fusion_plasmaEURprot.tsv
