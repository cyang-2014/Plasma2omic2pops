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
## update https://github.com/opain/Calculating-FUSION-TWAS-weights-pipeline/blob/master/OP_packaging_fusion_weights.R

# Create a file containing wgt.RDat file name, Gene ID, CHR, P0 and P1.
# Create file containing a list of the .wgt.RDat files
RDat_dir <- '/home/yangc/mendelianRandomization/TWAS/fusion/washu_files/test03_plasmaProtEURcisNtrans/s1_compute_weights/d02output/Output'
temp = list.files(path=RDat_dir, pattern="*.RDat")
temp_withPath<-paste('TWAS_Weights','/', temp, sep='')



system(paste0('mkdir /home/yangc/mendelianRandomization/TWAS/fusion/washu_files/test03_plasmaProtEURcisNtrans/s1_compute_weights/d02output/TWAS_weights_package'))
output_dir <- '/home/yangc/mendelianRandomization/TWAS/fusion/washu_files/test03_plasmaProtEURcisNtrans/s1_compute_weights/d02output/TWAS_weights_package'
output_name <- 'TWAS_Weights'
write.table(temp_withPath, paste(output_dir,'/',output_name,'.list',sep=''), col.names=F, row.names=F, quote=F)


# Create folder and insert .wgt.RDat files
system(paste('mkdir -p ',output_dir,'/',output_name,sep=''))
system(paste('cp ',RDat_dir,'/*.RDat ',output_dir,'/',output_name,'/',sep=''))


coordinate_file <- '/home/yangc/mendelianRandomization/TWAS/fusion/washu_files/test01_plasmaProtEUR/s1_compute_weights/Calculating-FUSION-TWAS-weights-pipeline/plasmaEUR_coordinate_file.txt'

pos_temp<-data.frame(   WGT=temp_withPath,
                        ID=gsub('.wgt.RDat','',gsub('Fetal_','',temp)))

Gene_coordinates_file<-read.table(coordinate_file, header=T, stringsAsFactors=F)
names(Gene_coordinates_file)<-c('CHR','start','end','ID')

pos_temp_2<-merge(pos_temp, Gene_coordinates_file[1:4], by='ID')
names(pos_temp_2)<-c('ID','WGT','CHR','P0','P1')
pos_temp_2$PANEL<-output_name
pos_temp_2<-pos_temp_2[c('PANEL','WGT','ID','CHR','P0','P1')]

pos_temp_2_sort<-pos_temp_2[order(pos_temp_2$CHR,pos_temp_2$P0),]
# pos_temp_2_sort


### maybe remove the covariate-files when calculating the weight from "s2"??


# Create a file containing the Gene ID, nsnps, hsq, hsq.se, hsq.pv, top1.r2, blup.r2, enet.r2, bslmm.r2, lasso.r2, top1.pv, blup.pv, enet.pv, bslmm.pv and lasso.pv (.profile)
profile_temp<-data.frame(       ID=gsub('.wgt.RDat','',temp),
                                nsnps=NA,
                                hsq=NA,
                                hsq.se=NA,
                                hsq.pv=NA,
                                top1.r2=NA,
                                blup.r2=NA,
                                enet.r2=NA,
                                top1.pv=NA,
                                blup.pv=NA,
                                enet.pv=NA)

pos_temp_2_sort$N<-NA
# pos_temp_2_sort

for(i in 1:nrow(pos_temp_2_sort)){
  # print(i)
  load(paste(RDat_dir,'/',pos_temp_2_sort$ID[i],'.wgt.RDat', sep=''))
  pos_temp_2_sort$N[[i]]<-N.tot
  profile_temp$nsnps[i]<-dim(snps)[1]
  profile_temp$hsq[i]<-hsq[1]
  profile_temp$hsq.se[i]<-hsq[2]
  profile_temp$hsq.pv[i]<-hsq.pv
  
  cv.performance<-data.frame(cv.performance)
  profile_temp$top1.r2[i]<-cv.performance$top1[1]
  profile_temp$blup.r2[i]<-cv.performance$blup[1]
  # profile_temp$bslmm.r2[i]<-cv.performance$bslmm[1]
  profile_temp$enet.r2[i]<-cv.performance$enet[1]
  # profile_temp$lasso.r2[i]<-cv.performance$lasso[1]
  
  profile_temp$top1.pv[i]<-cv.performance$top1[2]
  profile_temp$blup.pv[i]<-cv.performance$blup[2]
  # profile_temp$bslmm.pv[i]<-cv.performance$bslmm[2]
  profile_temp$enet.pv[i]<-cv.performance$enet[2]
  # profile_temp$lasso.pv[i]<-cv.performance$lasso[2]
}

nrow(pos_temp_2_sort) # 2285, not 812 loaded
nrow(profile_temp) # 2333


## note: not all proteins have a weight due to heritability estimation needs to be > 0 even without considering the h2_pvalue
system(paste0("ls /home/yangc/mendelianRandomization/TWAS/fusion/washu_files/test03_plasmaProtEURcisNtrans/s1_compute_weights/d02output/Output/*.wgt* | wc -l"))
## 2333
pos_temp_2_sort.dt <- data.table(pos_temp_2_sort)
profile_temp.dt <- data.table(profile_temp)

#### check the heritability for each weight calculation is significant or not
hist(profile_temp.dt$hsq.pv, breaks = 10)
profile_temp.dt[hsq.pv < 0.01]
profile_temp.dt[hsq.pv < 0.05]

profile_temp.dt[hsq.pv >= 0.5]
profile_temp.dt[hsq.pv >= 0.1]
## note: not all proteins have a non-NA weight per elastic-net! also, 48 proteins have missing weight per top1 & blup
profile_temp.dt[!is.na(enet.r2)]
profile_temp.dt[is.na(enet.r2)] # 1320 proteins with missing weight
table(profile_temp.dt[is.na(enet.r2)]$hsq.pv > 0.01)

profile_temp.dt[is.na(top1.r2)]
profile_temp.dt[is.na(blup.r2)]

write.table(pos_temp_2_sort, paste(output_dir,'/',output_name,'.pos',sep=''), col.names=T, row.names=F, quote=F)
write.table(profile_temp,paste(output_dir,'/',output_name,'.profile',sep=''), col.names=T, row.names=F, quote=F)

table(pos_temp_2_sort$CHR)



# Create a file comparing the different models
se <- function(x) sd(x)/sqrt(length(x))
BEST_r2<-c(	profile_temp$top1.r2[profile_temp$top1.r2 == max(profile_temp$top1.r2)],
            profile_temp$blup.r2[profile_temp$blup.r2 == max(profile_temp$blup.r2)],
            profile_temp$enet.r2[profile_temp$enet.r2 == max(profile_temp$enet.r2)])

profile_temp_r2<-profile_temp[c('top1.r2','blup.r2','enet.r2')]
for(k in 1:dim(profile_temp_r2)[1]){
  profile_temp_r2[k,][profile_temp_r2[k,] != max(profile_temp_r2[k,])]<-NA
}

sink(file = paste(output_dir,'/',output_name,'.profile.err',sep=''))
cat('Average hsq: ',mean(profile_temp$hsq),' ( ',se(profile_temp$hsq),' )


Average crossvalidation R2:
R2\tSE
top1\t',round(mean(profile_temp$top1.r2),3),'\t',round(se(profile_temp$top1.r2),5),'
blup\t',round(mean(profile_temp$blup.r2),3),'\t',round(se(profile_temp$blup.r2),5),'
enet\t',round(mean(profile_temp$enet.r2),3),'\t',round(se(profile_temp$enet.r2),5),'
BEST\t',round(max(BEST_r2),3),'
% Model is best:
top1:\t',round(sum(!is.na(profile_temp_r2$top1.r2))/length(profile_temp_r2$top1.r2)*100,1),'%
blup:\t',round(sum(!is.na(profile_temp_r2$blup.r2))/length(profile_temp_r2$blup.r2)*100,1),'%
enet:\t',round(sum(!is.na(profile_temp_r2$enet.r2))/length(profile_temp_r2$enet.r2)*100,1),'%
\n', sep='')
sink()



