library("VariantAnnotation")
library("vcfR")
library(snpStats)
library(R.utils)
library(data.table)
library(stringr)
library(dplyr)

###Get the variables we need
library(config)
Sys.setenv(R_CONFIG_ACTIVE = "300bcg")
config <- config::get()
for (x in names(config)){eval(parse(text=paste0(x,"='",config[[x]],"'")))}

source('processingFxnsGeneral.R')
source('plottingFxnsGeneral.R')
source('300bcg_generalFxns.R')

convertVcfToDosage <- function(folderWithVcfFiles,chrSel,fileSel,subFolderWithData){
  fixSampleNames <- function(oldNames){
    library(stringr)
    newNames = str_split_fixed(oldNames, '_',2)[,1]
    newNames = str_pad(string = newNames, width = 3, pad = '0', side = 'left')
    newNames = paste0('BCG',newNames)
    return(newNames)
  }
  
  print(fileSel)
  vcfFileName = file.path(folderWithVcfFiles,subFolderWithData,fileSel)
  
  #Extract info on the genotypes
  print(paste('saving info file for file',fileSel,sep=''))
  vcf2 <- read.vcfR(file = vcfFileName)
  vcfInfo = getFIX(vcf2, getINFO = TRUE)
  rm(vcf2)
  
  #The cases where no rs-id was assigned it will automatically add something like this in the dosages:
  #chr6:302973_C/G
  #chr6:445706_G/A
  rowsNoID = which(is.na(vcfInfo[,'ID']))
  if (length(rowsNoID)>0){
    vcfInfo[rowsNoID,'ID'] = paste0(vcfInfo[rowsNoID,'CHROM'],':',vcfInfo[rowsNoID,'POS'],'_',vcfInfo[rowsNoID,'REF'],'/',vcfInfo[rowsNoID,'ALT'])
  }
  
  #Now remove the few SNPs that were mapped to a different chromosome for some reason
  if (length(which(vcfInfo[,'CHROM']!=paste0('chr',chrSel)))>0){
    wrongChrom = vcfInfo[which(vcfInfo[,'CHROM'] != paste0('chr',chrSel)),,drop=F]
    idsToRemove = wrongChrom[,'ID']
    vcfInfo = vcfInfo[which(!(vcfInfo[,'ID'] %in% idsToRemove)),]
  }else{
    idsToRemove=NULL
  }
  
  #Now remap the positions to rs-ids for those that were updated
  fileWithNewRs = file.path(resourceDir,"dbSNP",paste0("human_9606_b151_GRCh38p7_SNPdb.filt_chr",chrSel,".tdf"))
  rsInfo.subset = fread(file = fileWithNewRs)
  rsInfo.subset$CHROM = paste0('chr',as.character(rsInfo.subset$CHROM))
  rsInfo.subset$POS = as.character(rsInfo.subset$POS)
  vcfInfo.merged = base::merge(vcfInfo, rsInfo.subset, by=c("CHROM", "POS","REF","ALT"), all.x=T, all.y=F)
  differentIds = which((vcfInfo.merged$ID.x!=vcfInfo.merged$ID.y) & !is.na(vcfInfo.merged$ID.y))
  vcfInfo.merged$ID.x = as.character(vcfInfo.merged$ID.x)
  vcfInfo.merged[differentIds,'ID.x'] = vcfInfo.merged[differentIds,'ID.y'] 
  vcfInfo.merged$ID.y <- NULL
  colnames(vcfInfo.merged) = gsub('\\.x$','',colnames(vcfInfo.merged))
  
  stopifnot(length(unique(vcfInfo[,'CHROM']))==1)
  
  #Save info in genotypes
  fileInfoName = paste('chr_',chrSel,'_300bcg_variantInfoAddedChrPos',sep='')
  fwrite(x = as.data.frame(vcfInfo), file = file.path(csvOutput,paste0(fileInfoName,'.txt')), sep = '\t', row.names = F, col.names = T)
  saveRDS(object =  as.data.frame(vcfInfo), file = file.path(rdsOutput,paste0(fileInfoName,'.RDS')))
  
  #Load the actual genotype data
  print(paste('load genotype data for chromosome',chrSel,sep=''))
  vcf <- readVcf(vcfFileName)
  
  #Keep only the ones with the correct chromosome
  toKeepIndices = which(as.vector(seqnames(rowRanges(vcf))==paste0('chr',chrSel)))
  vcf = vcf[toKeepIndices,]
  
  #(1) Hard calls
  mat <- genotypeToSnpMatrix(vcf, uncertain=F)
  hardCalls = as(mat$genotypes, 'numeric')
  hardCalls = t(hardCalls)
  oldNames = colnames(hardCalls)
  newNames = fixSampleNames(oldNames)
  colnames(hardCalls) = newNames
  hardCalls = hardCalls[,order(colnames(hardCalls))]
  
  #Check that the reference ajd alternative definitions are the same for both methods of loading the vcf data
  stopifnot(all.equal(as.character(mat$map$allele.1), vcfInfo[,'REF']))
  stopifnot(all.equal(as.character(unlist(mat$map$allele.2)), vcfInfo[,'ALT']))
  stopifnot(all.equal(as.character(unlist(mat$map$snp.names)), vcfInfo[,'ID']))
  
  #(2) Dosages
  dosages <- geno(vcf)$DS
  oldNames = colnames(dosages)
  newNames = fixSampleNames(oldNames)
  colnames(dosages) = newNames
  dosages = dosages[,order(colnames(dosages))]
  
  #Save the RDS and text files for both hard calls and dosages
  fwrite(x = as.data.frame(hardCalls), file = file.path(csvOutput,gsub('\\.vcf\\.gz','_hardCalls.csv',fileSel)), sep = '\t', row.names = T, col.names = T)
  saveRDS(object =  as.data.frame(hardCalls), file = file.path(rdsOutput,gsub('\\.vcf\\.gz','_hardCalls.RDS',fileSel)))
  
  fwrite(x = as.data.frame(dosages), file = file.path(csvOutput,gsub('\\.vcf\\.gz','_dosages.csv',fileSel)), sep = '\t', row.names = T, col.names = T)
  saveRDS(object =  as.data.frame(dosages), file = file.path(rdsOutput,gsub('\\.vcf\\.gz','_dosages.RDS',fileSel)))
  
  #Remove from mem
  rm(vcf)
  rm(mat)
  rm(dosages)
  rm(hardCalls)
}

#Inputs
subFolderWithData= 'mafHweOutlierFilteredGrCh38' #filtered data
folderWithVcfFiles= file.path(genetics_save_dir,"vcf") #vcf folder

#Save text files and RDS files
csvOutput=file.path(genetics_save_dir,'CSV',subFolderWithData)
rdsOutput=file.path(genetics_save_dir,'RDS',subFolderWithData)
mkdirs(csvOutput)
mkdirs(rdsOutput)

###The input argument should be for which chromosome to run this file
args = commandArgs(trailingOnly=TRUE)
chrSel = args[[1]]
#chrSel = 1 #Chromosome is set as an input argument --> run this script for each chromosone 1-22

fileSel = paste0('GRCh38_300BCG_chr_',chrSel,'_annotated_MAF0.1_HWE1e-5_noOutliers.recode.tags.vcf.gz')
convertVcfToDosage(folderWithVcfFiles,chrSel,fileSel,subFolderWithData)