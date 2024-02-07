#######################################################
## to create my own WEIGHTS for each protein-coding gene (6907 in total)
## using the pipeline built by Oliver Pain
# https://github.com/opain

## the github repo:
# https://github.com/opain/TWAS-plotter

rm(list = ls())
options(stringsAsFactors = FALSE)
library(reshape2)
library(dplyr)
library(data.table)
library(ggplot2)
library(cowplot)
library(Biobase)
library(gdata)
library(foreach)
library(ggrepel)
## use the plasma3170-proteomics data as a test


#################################################
## https://github.com/opain/TWAS-plotter/blob/master/TWAS-plotter.V1.0.r

# Create plotting function
TWAS_manhattan = function(dataframe, title=NULL, ylimit=max(abs(dataframe$TWAS.Z),na.rm=T)+1, Sig_Z_Thresh=qnorm(1-(0.05/length(dataframe$TWAS.Z))/2)) {
  
  d=dataframe[order(dataframe$CHR, dataframe$P0),]
  d=d[!is.na(d$TWAS.P),]
  
  d$pos=NA
  ticks=NULL
  lastbase=0
  numchroms=length(unique(d$CHR))
  
  for (i in unique(d$CHR)) {
    if (i==1) {
      d[d$CHR==i, ]$pos=d[d$CHR==i, ]$P0
    }	else {
      lastbase=lastbase+tail(subset(d,CHR==i-1)$P0, 1)
      d[d$CHR==i, ]$pos=d[d$CHR==i, ]$P0+lastbase
    }
    ticks=c(ticks, d[d$CHR==i, ]$pos[floor(length(d[d$CHR==i, ]$pos)/2)+1])
  }
  
  ticklim=c(min(d$pos),max(d$pos))
  
  mycols=rep(c("gray35","gray72"),60)
  
  d$Sig_Z_Thresh<-Sig_Z_Thresh
  
  d_sig<-d[which(abs(d$TWAS.Z) > d$Sig_Z_Thresh),]
  d_sig<-d_sig[rev(order(abs(d_sig$TWAS.Z))),]
  d_sig<-d_sig[!duplicated(d_sig$ID),]
  
  if(sum(d_sig$TWAS.Z > 0) > 0){
    d_sig_pos<-d_sig[d_sig$TWAS.Z > 0,]
  }
  
  if(sum(d_sig$TWAS.Z < 0) > 0){
    d_sig_neg<-d_sig[d_sig$TWAS.Z < 0,]
  }
  
  chr_labs<-as.character(unique(d$CHR))
  chr_labs[chr_labs == '19'| chr_labs == '21']<-' '
  
  if(dim(d_sig)[1] == 0){
    p<-ggplot(d,aes(x=pos,y=TWAS.Z,colour=factor(CHR))) +
      geom_point(size=0.5) +
      scale_x_continuous(name="Chromosome", breaks=ticks, labels=chr_labs) +
      scale_y_continuous(name='Z score',limits=c(-ylimit,ylimit)) +
      scale_colour_manual(values=mycols, guide=FALSE) +
      geom_hline(yintercept=0,colour="black") +
      geom_hline(yintercept=Sig_Z_Thresh,colour="blue") +
      geom_hline(yintercept=-Sig_Z_Thresh,colour="blue")
    
  } else {
    
    p<-ggplot(d,aes(x=pos,y=TWAS.Z,colour=factor(CHR))) +
      geom_point(size=0.5) +
      scale_x_continuous(name="Chromosome", breaks=ticks, labels=chr_labs) +
      scale_y_continuous(name='Z score',limits=c(-ylimit,ylimit)) +
      scale_colour_manual(values=mycols, guide=FALSE) +
      geom_hline(yintercept=0,colour="black") +
      geom_hline(yintercept=Sig_Z_Thresh,colour="blue") +
      geom_hline(yintercept=-Sig_Z_Thresh,colour="blue") +
      geom_point(data=d_sig, aes(x=pos,y=TWAS.Z), colour="red", fill='red', size=1.5)
    
    if(sum(d_sig$TWAS.Z > 0) > 0){
      p<-p+geom_text_repel(data=d_sig_pos, aes(x=pos,y=TWAS.Z, label=ID), colour='black', nudge_y=1, size=2.5, force=5, segment.alpha=0.25, ylim=c(Sig_Z_Thresh+0.1,NA))
    }
    
    if(sum(d_sig$TWAS.Z < 0) > 0){
      p<-p+geom_text_repel(data=d_sig_neg, aes(x=pos,y=TWAS.Z, label=ID), colour='black', nudge_y=-1, size=2.5, force=5, segment.alpha=0.25, ylim=c(NA,-Sig_Z_Thresh-0.1))
    }
  }
  
  p<-p + 	
    theme_cowplot() +
    theme(axis.text.x = element_text(angle=45, size=8, hjust=1))
  
  p
}

# Read in the TWAS data
## load in fusion output
DIRout <- '/home/yangc/mendelianRandomization/TWAS/fusion/washu_files/test03_plasmaProtEURcisNtrans/s2_perform_expression_imputation/eurADrisk2022/output/'
list.files(DIRout)

all22chr.twas <- foreach(i=1:22, .combine = 'rbind') %do% {
  
  singleCHR.twas <- fread(paste0(DIRout, 'eurADrisk2022fusion_plasmaEURprot.tsv.chr', i))
  return(singleCHR.twas)
}

chr6.MHC.twas <- fread(paste0(DIRout, 'eurADrisk2022fusion_plasmaEURprot.tsv.chr6.MHC'))
all22chr.twas.mhc <- rbindlist(list(all22chr.twas, chr6.MHC.twas))

coord.dt <- fread('/home/yangc/Rprojects/multiOmics/reports/forQTL_identification/PIGEON_proteomicsQC/geneLocations/afterQC_subsetSplitFeature6989human_annotations.csv')
nrow(coord.dt)
all22chr.twas.annot <- merge(all22chr.twas.mhc, coord.dt[, list(ID=Analytes, Target)], by = 'ID')
all22chr.twas.annot[, ID := Target]
all22chr.twas.annot[is.na(TWAS.Z)]
all22chr.twas.annot[is.na(TWAS.Z)][NSNP > 0]

twas<-data.frame(all22chr.twas.annot)

TWAS_manhattan(dataframe=twas, Sig_Z_Thresh=3)

