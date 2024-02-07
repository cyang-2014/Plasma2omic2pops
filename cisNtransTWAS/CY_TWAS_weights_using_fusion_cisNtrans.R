#!/usr/bin/Rscript
library("optparse")

option_list = list(
  make_option("--PLINK_prefix", action="store", default=NA, type='character',
              help="Prefix of PLINK files for the target sample [required]"),
  make_option("--PLINK_prefix_chr", action="store", default=NA, type='character',
              help="Prefix of per chromosome PLINK files for the target sample [required]"),
  make_option("--phenotype_file", action="store", default=NA, type='character',
              help="File name for normalised and adjusted expression data [required]"),
  make_option("--gene_name", action="store", default=NA, type='character',
              help="Name of gene or transcript to be processed [required]"),
  make_option("--plink", action="store", default=NA, type='character',
              help="Path to PLINK [required]"),
  make_option("--gcta", action="store", default=NA, type='character',
              help="Path to gcta_nr_robust binary [required]"),
  make_option("--gemma", action="store", default=NA, type='character',
              help="Path to gemma binary [required]"),
  make_option("--ld_ref_dir", action="store", default=NA, type='character',
              help="FUSION LD reference directory [required]"),
  make_option("--fusion_software", action="store", default=NA, type='character',
              help="FUSION software directory [required]"),
  make_option("--output_dir", action="store", default=NA, type='character',
              help="Directory name for the output [required]"),
  make_option("--extractSNP_dir", action="store", default=NA, type='character',
              help="Directory name for the extracted variants in cis and trans per the same feature [required]")
)

opt = parse_args(OptionParser(option_list=option_list))
system(paste('mkdir ',opt$output_dir,'/Output',sep=''))

root_dir<-getwd()

library(data.table)


# Extract gene from phenotype file
system(paste('mkdir ',opt$output_dir,'/temp',sep=''))
system(paste('awk -f /home/yangc/mendelianRandomization/TWAS/fusion/washu_files/test01_plasmaProtEUR/s1_compute_weights/Calculating-FUSION-TWAS-weights-pipeline/d01scripts/t.awk c1=FID c2=IID c3=',opt$gene_name,' ',opt$phenotype_file,' > ',opt$output_dir,'/temp/temp_',opt$gene_name,'.pheno',sep=''))


# Using PLINK, extract variants for the cis+trans region of each feature
err_1<-system(paste(opt$plink,' --bfile ',opt$PLINK_prefix,' --make-bed --pheno ',opt$output_dir,'/temp/temp_',opt$gene_name,'.pheno',
                    ' --out ', opt$output_dir,'/temp/temp_',opt$gene_name,
                    ' --geno 0.02 ', '--extract ',
                    opt$extractSNP_dir, '/', opt$gene_name, '_EUR.txt', sep='')) # note: this "_EUR" string also needs to be changed per analysis


if (err_1 == 13) {
write.table('No SNPs within the cis and trans region given this feature', paste(opt$output_dir,'/Output/',opt$gene_name,'.err',sep=''), col.names=F, row.names=F, quote=F)
} else {

# Using FUSION, calculate the weights for the genes expression using subset of genotypic data.
setwd(paste(opt$output_dir,'/temp', sep=''))
system(paste('ln -s ./ output', sep=''))

#system(paste('Rscript ',opt$fusion_software,'/FUSION.compute_weights.R --bfile ', opt$output_dir,'/temp/temp_',opt$gene_name,' --tmp temp_',opt$gene_name,'.tmp --out ', opt$output_dir,'/Output/',opt$gene_name,' --verbose 2 --save_hsq --PATH_gcta ',opt$gcta,' --PATH_gemma ',opt$gemma,' --PATH_plink ',opt$plink,' --models top1,lasso,enet,blup', sep=''))

system(paste('Rscript ',opt$fusion_software,'/FUSION.compute_weights.R --bfile ', opt$output_dir,'/temp/temp_',opt$gene_name,' --tmp temp_',opt$gene_name,'.tmp --out ', opt$output_dir,'/Output/',opt$gene_name,' --verbose 2 --hsq_p 1 --save_hsq --PATH_gcta ',opt$gcta,' --PATH_gemma ',opt$gemma,' --PATH_plink ',opt$plink,' --models top1,enet,blup', ' --crossval 2', sep=''))

}
# Delete the temporary files
#system(paste('rm ',opt$output_dir,'/temp/temp_',opt$gene_name,'*', sep=''))
