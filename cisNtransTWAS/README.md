# FUSION/TWAS with cis and trans regions of molecular traits


## fusion/TWAS was developed by Sasha Gusev et al at https://github.com/gusevlab/fusion_twas
## the pipeline is also inspired by Oliver Pain at https://github.com/opain/Calculating-FUSION-TWAS-weights-pipeline
## the pipeline only included the cis-regions of the gene; and this pipeline was considering both cis and trans regions given the same gene/protein.
## this cis+trans design futher made the metabolite feasible, as no metabQTL was strictly defined as cis metabQTL, compared to eQTL/pQTL.
## Chengran Yang


#########################
# partA Create weights
#########################
### the core script "CY_TWAS_weights_using_fusion_cisNtrans.R"
#### arguments:
#### --PLINK_prefix: the prefix of the LD reference file in plink format
#### --phenotype_file: feature by sample matrix (with FID, IID as the first two columns, the rest columns as the featureID, rows are sampleID), similar input format for plink2 glm association test, but better regress out the covariate-matrix
#### --gene_name: featureID, not necessarily the gene-name, can be protein-name, metabolite-name, etc
#### --plink: path of the plink1.9 software
#### --gcta: path of the GCTA package
#### --gemma: path of the gemma package 
#### --ld_ref_dir: LD reference file in plink format
#### --fusion_software: local path after downloading the fusion package 
#### --output_dir: output directory path
#### --extractSNP_dir: storing both the cis and trans regions of the study-wide significant variants for each feature, unique input compared to cis only framework 


### step1, prepare the input 
## in Rstudio, update the script "s01_fusion_twas_prep_input_files_option2_EUR.R"


### step-2, loop all features after assigning all the customized arguments
Rscript /home/yangc/Rprojects/multiOmics/code/phase14_TWASmethods/PWAS_GusevLab/option2_cisPlusTrans/p01_EURprot/s2_fusion_twas_o2_process_allProteins_EUR.R &


### step-3, combine all weights in one folder
## in Rstudio, process "s03_fusion_twas_o2_combineWeights_fromAllProteins_EUR.R" after updating the path for the files to be read


#########################
# partB Perform association tests
#########################
## I updated the "FUSION.assoc_test.R" to "CY.FUSION.assoc_test.R"
# the idea is to break the single-chromosome requirement per protein-disease association

## process CY.FUSION.assoc_test.R using the script "s04_new_assocTest_noCHR.bash"
bash s04_new_assocTest_noCHR.bash

#########################
# partC Visualize the TWAS results with Miami plot
#########################
## in Rstudio, update the script "s05_fusion_twas_o2_assocTest_fromAllProteins_EUR_eurADrisk2022_Miamiplot.R"
