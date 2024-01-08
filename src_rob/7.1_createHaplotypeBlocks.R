library("R.utils")
library('openxlsx')
library(ggplot2)
library(stringr)
library(data.table)
library(dplyr)
library(stats)
library(gpart)

print('starting')

###Get the variables we need
library(config)
Sys.setenv(R_CONFIG_ACTIVE = "300bcg")
config <- config::get()
for (x in names(config)){eval(parse(text=paste0(x,"='",config[[x]],"'")))}

source('processingFxnsGeneral.R')
source('plottingFxnsGeneral.R')
source('300bcg_generalFxns.R')

outputFolderForHaplotype = file.path(data_dir_300BCG,'genetics','haplotypeBlocks_BigLD')
mkdirs(outputFolderForHaplotype)

#Define haplotype blocks
#load genetics data
subFolderWithData = "mafHweOutlierFilteredGrCh38"
rdsOutput=file.path(genetics_save_dir,'RDS',subFolderWithData)

#The input argument is the chromsome, so it should be run for chromosome 1 to 22
args = base::commandArgs(trailingOnly=TRUE)
chrSel = args[[1]]

saveFileGeneticsChrom = file.path(rdsOutput,paste0("GRCh38_300BCG_chr_",chrSel,"_annotated_MAF0.1_HWE1e-5_noOutliers.recode.tags_hardCalls.RDS"))
geno = readRDS(saveFileGeneticsChrom)
geno = t(geno)

variantInfo = readRDS(file.path(rdsOutput,paste0("chr_",chrSel,"_300bcg_variantInfoAddedChrPos.RDS")))
variantInfo = variantInfo[which(variantInfo$CHROM==paste0("chr",chrSel)),]

geno = geno[,as.character(variantInfo$ID)]
##   chrN        rsID       bp
## 1   21  rs59504523 16040232
## 2   21 rs193065234 16040492
## 3   21 rs117862161 16040674
## 4   21 rs185637933 16040851

SNPinfo = variantInfo[,c('CHROM','ID','POS')]
colnames(SNPinfo) = c("chrN","rsID","bp")
SNPinfo$bp = as.numeric(as.character(SNPinfo$bp))
SNPinfo$chrN = as.numeric(gsub('^chr','',as.character(SNPinfo$chrN)))

res1 = gpart::BigLD(geno = geno, SNPinfo = SNPinfo) 
fwrite(res1, file = file.path(outputFolderForHaplotype,paste0("GRCh38_300BCG_chr_",chrSel,"_annotated_MAF0.1_HWE1e-5_noOutliers.recode.tags.haplo")))