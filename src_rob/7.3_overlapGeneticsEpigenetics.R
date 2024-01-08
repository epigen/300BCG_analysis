library("R.utils")
library('openxlsx')
library(ggplot2)
library(stringr)
library(data.table)
library(dplyr)
library(stats)

# Define a few helper functions
#https://jef.works/blog/2016/12/06/mapping-snps-and-peaks-to-genes-in-R/
#' Convert from string to range
#' 
#' @param pos A vector of strings ex. chr1 2938302 2938329
#' @param delim Delimiter for string splitting
#' @param region Boolean of whether region or just one position
#'
#' @returns Dataframe of ranges
#' 
string2range <- function(pos, delim=' ', region=TRUE) {
  posp <- as.data.frame(do.call(rbind, strsplit(pos, delim)))
  posp[,1] <- posp[,1]
  posp[,2] <- as.numeric(as.character(posp[,2]))
  if(region) {
    posp[,3] <- as.numeric(as.character(posp[,3]))
  } else {
    posp[,3] <- posp[,2]
  }
  return(posp)
}

#' Convert from ranges to GRanges
#' 
#' @param df Dataframe with columns as sequence name, start, and end
#' 
#' @returns GRanges version 
#' 
range2GRanges <- function(df) {
  require(GenomicRanges)
  require(IRanges)
  gr <- GenomicRanges::GRanges(
    seqnames = df[,1],
    ranges=IRanges(start = df[,2], end = df[,3])
  )
  return(gr)
}

print('starting')

###Get the variables we need
library(config)
Sys.setenv(R_CONFIG_ACTIVE = "300bcg")
config <- config::get()
for (x in names(config)){eval(parse(text=paste0(x,"='",config[[x]],"'")))}

source('300bcg_generalFxns.R')

peakInfoFile.haplo = file.path(folderPeaksExpanded,'haplotypeExpandedAtacSeqRegions_total.csv')

#Take the best SNP value (i.e. lowest pvalue) for each region per cytokine-stim combination
#Row-bind them all together
#After that do fisher exact test

#Examples of how to run the script
# sbatch 300bcg_overlapGeneticsEpigenetics.batch cytoData.good.log2 V1
# sbatch 300bcg_overlapGeneticsEpigenetics.batch circMarkerData.log2 V1
# sbatch 300bcg_overlapGeneticsEpigenetics.batch cytoData.good.fc V3
# sbatch 300bcg_overlapGeneticsEpigenetics.batch cytoData.good.fc V2
# Rscript 300bcg_overlapGeneticsEpigenetics.R cytoData.good.log2 V1
# Rscript 300bcg_overlapGeneticsEpigenetics.R circMarkerData.log2 V1
# Rscript 300bcg_overlapGeneticsEpigenetics.R cytoData.good.fc V3
###Example arguments
#args = c("cytoData.good.fc",'V3')
#args = c("cytoData.good.log2",'V1')
#args =  c('circMarkerData.log2','V1')

#Load the regions/peaks
####Here select either the 'old' peakinfo or the new expanded one based on haplotype blocks
peakInfo = fread(peakInfoFile.haplo)
peakInfo.orig = fread(peakInfoFile)
#Load all the data and metadata
allRelevantDataList = loadLukasMetaDataPbmc(allMetaDataPath)
list2env(allRelevantDataList, .GlobalEnv)

args = base::commandArgs(trailingOnly=TRUE)
#args = c("cytoData.good.log2",'V1')
dataNameSel = args[[1]]
timePointSel = args[[2]]

#Extract the data of interest
dataSel = allRelevantDataList[[dataNameSel]]

#Create a list of folders where to load ATAC data
atacResultsFolderList = list()
atacResultsFolderList[["cytoData.good.log2"]] = file.path(folderWithBaselineAtacResults,"cytokines")
atacResultsFolderList[["cytoData.good.fc"]] = file.path(folderWithTrainingAtacResults,"cytokines")
atacResultsFolderList[["circMarkerData.log2"]] = file.path(folderWithBaselineAtacResults,"CM")
atacResultsFolderList[["circMarkerData.fc"]] = file.path(folderWithTrainingAtacResults,"CM")

#Get the folder with the genetics results
subqtlDirOverall = list.files(include.dirs = T, path = qtlDirOverall, pattern = paste0(dataNameSel,'__',timePointSel))
fullFolderGenetics = file.path(qtlDirOverall,subqtlDirOverall,"perChrom")

#Folder to save overlap results
folderToSaveComparisonFile = file.path(data_dir_300BCG,'genetics','atacGeneticsCompared')
mkdirs(folderToSaveComparisonFile)

redoOverlap=T
if (redoOverlap){
  totalMergedData = NULL
  ##Iterate over all cytokines/variables
  for (variableSel in colnames(dataSel)){
    print(variableSel)
    variableSelShort = gsub('.fc$','',variableSel)
    if (variableSelShort %in% cytosToExclude){next}
    #Load the epigenetics results
    #grep the file of interest
    folderWithAtacResults = atacResultsFolderList[[dataNameSel]]
    if (endsWith(dataNameSel,'.fc')){
      #V1 = d0; V2 = d14; V3 = d90
      selFile = list.files(folderWithAtacResults, pattern = paste0('.*',timePointSel,'.*',gsub('\\.fc$','',variableSel),'.*\\.csv\\.gz$'))
    }else{
      selFile = list.files(folderWithAtacResults, pattern = paste0('.*',variableSel,'.*\\.csv\\.gz$'))
    }
    atacResults = fread(file.path(folderWithAtacResults,selFile))
    #add peak information
    atacResults.complete = base::merge(atacResults, peakInfo, by.x = 'V1', by.y = 'peak_id')
    
    #Create Granges object ATAC
    library(GenomicRanges)
    peaks.granges <- GRanges(
      seqnames = atacResults.complete$chr,
      ranges = IRanges(start = atacResults.complete$start, end = atacResults.complete$end, names = atacResults.complete$V1),
      score = atacResults.complete$F.p.value)

    #Load the genetics results
    #Grep the filenames
    #Go over all chromosomes
    allQtlFiles = list.files(path = fullFolderGenetics, pattern = paste0('^',variableSel,"_chr.*_withInfo.cqtls"))
    qtlResPerVar = NULL
    for (qtlFileSel in allQtlFiles){
      print(qtlFileSel)
      qtlResChrom = fread(file = file.path(fullFolderGenetics,qtlFileSel))
      if (is.null(qtlResPerVar)){
        qtlResPerVar = qtlResChrom
      }else{
        qtlResPerVar = rbind(qtlResPerVar, qtlResChrom)
      }
    }
    stopifnot(sum(is.na(qtlResPerVar$snps))==0)
    stopifnot(sum(is.na(qtlResPerVar$CHROM))==0) 
    
    #Create GRanges object genetics
    snps.granges <- GRanges(
      seqnames = qtlResPerVar$CHROM,
      ranges = IRanges(start = qtlResPerVar$POS, width = 1, names = qtlResPerVar$snps),
      score = qtlResPerVar$pvalue)
    
    #Calculate overlap between genetics and epigenetics
    overlap <- GenomicRanges::findOverlaps(peaks.granges,snps.granges)
    
    #Go over each region and select the best/lowest pvalue from the cQTL data
    overlap.dt = as.data.table(overlap)
    allRegionsIndicesWithSnp = unique(overlap.dt$queryHits)
    setkey(overlap.dt, queryHits)
    
    #matchedSnps.dt = NULL
    selList = list()
    for (regionIndex in allRegionsIndicesWithSnp){
      if(regionIndex %% 1000==0){print(regionIndex)}
      regionSel = atacResults.complete[[regionIndex,'V1']]
      overlap.sel = overlap.dt[J(regionIndex),]
      indicesOfSnps = overlap.sel$subjectHits
      snpSubset = qtlResPerVar[indicesOfSnps,]
      snpSubset = snpSubset[which.min(snpSubset$pvalue),]
      snpSubset$regionMatched = regionSel
      
      selList[[regionIndex]] = snpSubset
    }
    folderToSaveCheck = file.path(folderToSaveComparisonFile,'someChecks')
    mkdirs(folderToSaveCheck)
    saveRDS(selList, file = file.path(folderToSaveCheck,paste0('listOfRegions_',dataNameSel,'_geneticsAtacCompare.RDS')))
    matchedSnps.dt = rbindlist(selList)
    
    #Now bind snps and ATAC regions results together
    matchedSnps.dt.forMerge = matchedSnps.dt[,c("snps","gene","statistic","beta","pvalue","regionMatched","CHROM","POS","REF","ALT")]
    colnames(matchedSnps.dt.forMerge) = paste0(colnames(matchedSnps.dt.forMerge),'_qtl')
    colsToSel = c('V1',"chr","start","end",
                  colnames(atacResults.complete)[grep('(t|p\\.value)\\..*', colnames(atacResults.complete))])
    atacResults.complete.formerge = atacResults.complete[,..colsToSel]
    colnames(atacResults.complete.formerge) = gsub('\\.(LFC|CYTO|CM).*','',colnames(atacResults.complete.formerge))
    colnames(atacResults.complete.formerge) = paste0(colnames(atacResults.complete.formerge),'_atac')
    
    mergedResults = base::merge(matchedSnps.dt.forMerge,atacResults.complete.formerge,by.x='regionMatched_qtl',by.y='V1_atac')
    saveRDS(mergedResults,file = file.path(folderToSaveComparisonFile, paste0(dataNameSel,'_',variableSel,'_geneticsAtacCompare.RDS')))
    
    #Combine into a full df with all the data (basically combine all cytokines in a single df)
    if (is.null(totalMergedData)){
      totalMergedData = mergedResults
    }else{
      totalMergedData = rbind(mergedResults,totalMergedData)
    }
  }
  saveRDS(totalMergedData,file = file.path(folderToSaveComparisonFile, paste0('totalDf_',dataNameSel,'_geneticsAtacCompare.RDS')))
}else{
  totalMergedData = readRDS( file.path(folderToSaveComparisonFile, paste0('totalDf_',dataNameSel,'_geneticsAtacCompare.RDS')))
}

#Create variable that counts both the region and cytokine
totalMergedData$cytoRegionCombi = paste(totalMergedData$regionMatched_qtl, totalMergedData$gene_qtl, sep='_')

#Define the groups we want to look at
cytokineGroups = defineSubgroupsCytokines()

if (dataNameSel %in% c("cytoData.good.fc")){
  subGroupsToTest = c("all_cyto_noLact")
}else if (dataNameSel == "cytoData.good.log2"){
  subGroupsToTest = c("all_cyto_noLact")
}else if (dataNameSel=="circMarkerData.log2"){
  subGroupsToTest = c("all")
}

# create style, in this case bold header
header_st <- openxlsx::createStyle(textDecoration = "Bold")  

for (subGroupSel in subGroupsToTest){
  #Take the subgroup we are interested in
  print(subGroupSel)
  cytosSel = cytokineGroups[[dataNameSel]][[subGroupSel]]
  
  totalMergedData.sel = totalMergedData[which(totalMergedData$gene_qtl %in% cytosSel),]
  
  #Check overlap for a marginally significant threshold
  pValCutoff = 1e-3
  signifHitsAtac = totalMergedData.sel$cytoRegionCombi[totalMergedData.sel$p.value_atac<=pValCutoff]
  signifHitsSnp = totalMergedData.sel$cytoRegionCombi[totalMergedData.sel$pvalue_qtl<=pValCutoff]
  nonSignifAtac = totalMergedData.sel$cytoRegionCombi[totalMergedData.sel$p.value_atac>pValCutoff]
  nonSignifSnp = totalMergedData.sel$cytoRegionCombi[totalMergedData.sel$pvalue_qtl>pValCutoff]
  
  signifInBoth = length(intersect(signifHitsSnp, signifHitsAtac))
  signifInJustSnp = length(intersect(signifHitsSnp, nonSignifAtac))
  signifInJustAtac = length(intersect(signifHitsAtac, nonSignifSnp))
  signifInNeither = length(intersect(nonSignifAtac, nonSignifSnp))
  #Create matrix for fisher test
  matrixSel = matrix(c(signifInBoth, signifInJustAtac,  signifInJustSnp, signifInNeither), nrow = 2, dimnames =
                       list(c("signifSnp", "nonSignifSnp"),
                            c("signifInAtac", "nonSignifAtac")))
  
  #Do fisher test
  result = fisher.test(matrixSel, alternative="greater")
  
  #Save the overlap
  totalMergedData.signifInBoth = totalMergedData.sel[which(totalMergedData.sel$cytoRegionCombi %in% intersect(signifHitsSnp, signifHitsAtac)),]
  totalMergedData.signifInBoth = base::merge(totalMergedData.signifInBoth,peakInfo,by.x='regionMatched_qtl',by.y='peak_id',all.x=T)
  
  #For now I add the single p-value to all rows. However, it is a single p-value for the total dataframe, so highly redundant, but just a way to save it in the same Excel file.
  totalMergedData.signifInBoth$totalOverlapPval = result$p.value
  totalMergedData.signifInBoth$percentageOfDarsSignifSnps = 100*nrow(totalMergedData.signifInBoth)/sum(totalMergedData$p.value_atac<=1e-3)
  
  #add original peak info as a names peakInfo.orig
  peakInfo.orig[[definitionsForSupplement[['region']]]] = paste(peakInfo.orig$chr, peakInfo.orig$start, peakInfo.orig$end, sep= '_')
  selCols = c('peak_id',definitionsForSupplement[['region']])
  peakInfo.orig.sel = peakInfo.orig[,..selCols]
  totalMergedData.signifInBoth = base::merge(totalMergedData.signifInBoth,peakInfo.orig.sel, 
                                            by.x = "regionMatched_qtl", by.y = "peak_id")
  
  ##Save both the raw file and the actual supplementary file
  ##Save this as a raw table
  folderForRawTable = file.path(suppTableDir,'xlsx_large','overlapGenEpigen')
  mkdirs(folderForRawTable)
  fileNameXlsx = file.path(folderForRawTable,paste0('TableS7_',tableAdditionToFileNames[['TableS7']],'_pt2.xlsx'))
  if (!(file.exists(fileNameXlsx))){
    wb <- openxlsx::createWorkbook()
  }else{
    wb <- openxlsx::loadWorkbook(fileNameXlsx)
  }

  sheetName = paste('overlap',dataTypeNameMapShort[[dataNameSel]],mapVisitShortNameToFinal[[timePointSel]],sep='_')

  if (!(sheetName %in% wb$sheet_names)){
    openxlsx::addWorksheet(wb = wb, sheet = sheetName)
  }

  openxlsx::writeData(wb = wb, sheet = sheetName, x = totalMergedData.signifInBoth, colNames = T, rowNames = F, headerStyle = header_st)
  openxlsx::saveWorkbook(wb = wb, file = fileNameXlsx, overwrite = T)
  
  ##Save this as a supplemental table
  ##Take out the most important info and map
  totalMergedData.signifInBoth = as.data.frame(totalMergedData.signifInBoth)
  selDataForSupp = totalMergedData.signifInBoth[,c('gene_qtl','snps_qtl','CHROM_qtl', 'POS_qtl', 'REF_qtl', 'ALT_qtl','statistic_qtl','pvalue_qtl',
                                                   definitionsForSupplement[['region']] ,'chr','start','end','t_atac','p.value_atac')]
  colnames(selDataForSupp) = c(definitionsForSupplement[[paste(dataNameSel,'LONG',sep='_')]],
                               paste0('Identifier',' (', definitionsForSupplement[['snp']],')'),
                               paste0(definitionsForSupplement[['chr']],' (', definitionsForSupplement[['snp']],')'),
                               paste0(definitionsForSupplement[['pos']],' (', definitionsForSupplement[['snp']],')'),
                               paste0(definitionsForSupplement[['ref']],' (', definitionsForSupplement[['snp']],')'),
                               paste0(definitionsForSupplement[['alt']],' (', definitionsForSupplement[['snp']],')'),
                               paste0(definitionsForSupplement[['statistic']],' (', definitionsForSupplement[['snp']],')'),
                               paste0(definitionsForSupplement[['pval']],' (', definitionsForSupplement[['snp']],')'),
                               paste0(definitionsForSupplement[['region']],' (', definitionsForSupplement[['ATAC']],')'), #chr1_1_2
                               paste0(definitionsForSupplement[["chr_haplotype"]],' (', definitionsForSupplement[['ATAC']],')'),
                               paste0(definitionsForSupplement[["start_haplotype"]],' (', definitionsForSupplement[['ATAC']],')'),
                               paste0(definitionsForSupplement[["end_haplotype"]],' (', definitionsForSupplement[['ATAC']],')'),
                               paste0(definitionsForSupplement[['statistic']],' (', definitionsForSupplement[['ATAC']],')'),
                               paste0(definitionsForSupplement[['pval']],' (', definitionsForSupplement[['ATAC']],')'))
  
  stopifnot(sum(as.numeric(str_split_fixed(selDataForSupp[[paste0(definitionsForSupplement[['region']],' (', definitionsForSupplement[['ATAC']],')')]],'_',3)[,2])>=
                  selDataForSupp[,paste0(definitionsForSupplement[["start_haplotype"]],' (', definitionsForSupplement[['ATAC']],')')])==nrow(selDataForSupp))
  stopifnot(sum(as.numeric(str_split_fixed(selDataForSupp[[paste0(definitionsForSupplement[['region']],' (', definitionsForSupplement[['ATAC']],')')]],'_',3)[,3])<=
                  selDataForSupp[,paste0(definitionsForSupplement[["end_haplotype"]],' (', definitionsForSupplement[['ATAC']],')')])==nrow(selDataForSupp))
  
  #map cytokines
  selDataForSupp[[definitionsForSupplement[[paste(dataNameSel,'LONG',sep='_')]]]] = gsub('\\.fc$','',selDataForSupp[[definitionsForSupplement[[paste(dataNameSel,'LONG',sep='_')]]]])
  selDataForSupp[[definitionsForSupplement[[paste(dataNameSel,'LONG',sep='_')]]]] = mapCytoNamesForSupplement(selDataForSupp[[definitionsForSupplement[[paste(dataNameSel,'LONG',sep='_')]]]], dataNameSel)

  folderToSaveXlsx = file.path(suppTableDir,'xlsx_small','raw')
  mkdirs(folderToSaveXlsx)
  fileNameXlsx = file.path(folderToSaveXlsx,paste0('TableS7_',tableAdditionToFileNames[['TableS7']],'_pt2.xlsx'))
  if (!(file.exists(fileNameXlsx))){
    wb <- openxlsx::createWorkbook()
  }else{
    wb <- openxlsx::loadWorkbook(fileNameXlsx)
  }
  
  sheetName = paste('overlapSnpAtac',definitionsForSupplement[[dataNameSel]],definitionsForSupplement[[timePointSel]],sep='_')
  
  if (!(sheetName %in% wb$sheet_names)){
    openxlsx::addWorksheet(wb = wb, sheet = sheetName)
  }
  
  openxlsx::writeData(wb = wb, sheet = sheetName, x = selDataForSupp, colNames = T, rowNames = F, headerStyle = header_st)
  openxlsx::saveWorkbook(wb = wb, file = fileNameXlsx, overwrite = T)
}
