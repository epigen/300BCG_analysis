library("R.utils")
library('openxlsx')
library(ggplot2)
library(stringr)
library(data.table)
library(dplyr)
library(stats)

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

#Load the regions/peaks
peakInfo = fread(peakInfoFile)

##Redefine the regions of the ATAC based on any haplotypes they overlap with
library(GenomicRanges)
peaks.granges.total = NULL
for (chrSel in rev(c(1:22))){
  print(chrSel)
  #Create Granges object ATAC for this chromosome
  peakInfoSel = peakInfo[which(peakInfo$chr == paste0('chr',chrSel)),]
  peaks.granges <- GRanges(
    seqnames = factor(peakInfoSel$chr),
    ranges = IRanges(start = peakInfoSel$start, end = peakInfoSel$end, names = peakInfoSel$peak_id))
  
  peaks.granges.df = as.data.frame(peaks.granges)
  
  #Create Granges for haplotype regions
  fileSel = file.path(outputFolderForHaplotype,paste0("GRCh38_300BCG_chr_",chrSel,"_annotated_MAF0.1_HWE1e-5_noOutliers.recode.tags.haplo"))
  haploTypeRegions = fread(fileSel)
  haploTypeRegions.granges <- GRanges(
    seqnames = factor(paste0('chr',haploTypeRegions$chr)),
    ranges = IRanges(start = haploTypeRegions$start.bp, end = haploTypeRegions$end.bp, names = haploTypeRegions$start.index))
  haploTypeRegions.granges.df = as.data.frame(haploTypeRegions.granges)
  
  overlappingRegions = findOverlaps(peaks.granges, haploTypeRegions.granges) 
  overlappingRegions.df = as.data.frame(overlappingRegions)
  
  #Now iterate over all atac-seq regions and expand them if they overlap a haplotype region
  peaks.granges.df.expanded = peaks.granges.df
  for (indexAtac in unique(overlappingRegions.df[,'queryHits'])){
    indicesHaplo = overlappingRegions.df[which(overlappingRegions.df[,'queryHits']==indexAtac),'subjectHits']
    
    startAtac = peaks.granges.df[indexAtac,'start']
    endAtac = peaks.granges.df[indexAtac,'end']
    
    startHaplo = min(haploTypeRegions.granges.df[indicesHaplo,'start'])
    endHaplo = max(haploTypeRegions.granges.df[indicesHaplo,'end'])
    
    totalStart = min(startAtac, startHaplo)
    totalEnd = max(endAtac,endHaplo)
    
    peaks.granges.df.expanded[indexAtac,'start'] = totalStart
    peaks.granges.df.expanded[indexAtac,'end'] = totalEnd
  }
  
  peaks.granges.df.expanded$chr = peaks.granges.df.expanded$seqnames
  peaks.granges.df.expanded$length = peaks.granges.df.expanded$end-peaks.granges.df.expanded$start
  peaks.granges.df.expanded$peak_id = rownames(peaks.granges.df.expanded)
  peaks.granges.df.expanded.sel = peaks.granges.df.expanded[,c('peak_id','chr','start','end')]  
  

  mkdirs(folderPeaksExpanded)
  fwrite(x = peaks.granges.df.expanded.sel, row.names = F , file = file.path(folderPeaksExpanded,paste0('haplotypeExpandedAtacSeqRegions_chr',chrSel,'.csv')))
  
  if (is.null(peaks.granges.total)){
    peaks.granges.total = peaks.granges.df.expanded.sel
  }else{
    peaks.granges.total = rbind.data.frame(peaks.granges.total, peaks.granges.df.expanded.sel)
  }
  
}
fwrite(x = peaks.granges.total, row.names = F , file = file.path(folderPeaksExpanded,paste0('haplotypeExpandedAtacSeqRegions_total.csv')))