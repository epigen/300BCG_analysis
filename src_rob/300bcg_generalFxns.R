#Define some cytokines we do not want to include in the analysis
cytosToExclude = c("C.albicans.yeast_24h_PBMC_IFNg","C.albicans.yeast_24h_PBMC_lactate","MTB_24h_PBMC_lactate")

##Load fonts that might not be present (optional)
#place these in a folder

# if (Sys.info()['sysname']!="Darwin"){
#   library(extrafont)
#   #Helvetica font
#   #This should be present in the specified folder below
#   font_import(paths = '/home/rterHorst/Fonts', prompt=F)
# }

library(Cairo)
CairoFonts(regular="Helvetica:style=Regular",
           bold="Helvetica:style=Bold",
           italic="Helvetica:style=Italic",
           bolditalic="Helvetica:style=Bold Italic,BoldItalic")

mapVisitShortNameToFinal = list()
mapVisitShortNameToFinal[['V1']] = 'd0'
mapVisitShortNameToFinal[['V2']] = 'd14'
mapVisitShortNameToFinal[['V3']] = 'd90'

distanceNameMapShort = list()
distanceNameMapShort[['GENE_AND_DISTAL_10kb']] = 'gene'
distanceNameMapShort[['TSS_PROXIMAL']] = 'tss'

pathwaySetNameMapShort = list()
pathwaySetNameMapShort[["KEGG_2019_Human_min15_max500"]] = 'Kegg'
pathwaySetNameMapShort[["GO_Biological_Process_2018_min15_max500"]] = 'Go'

dataTypeNameMapShort = list()
dataTypeNameMapShort[['cytoData.good.log2']] = 'cytokine'
dataTypeNameMapShort[['circMarkerData.log2']] = 'inflaMark'
dataTypeNameMapShort[['cytoData.good.fc']] = 'cytokineFC'
dataTypeNameMapShort[['circMarkerData.fc']] = 'inflaMarkFC'
dataTypeNameMapShort[['cytoData.good.fcCorr']] = 'cytokineFC.cor'
dataTypeNameMapShort[['circMarkerData.fcCorr']] = 'inflaMarkFC.cor'

definitionsForSupplement = list()
definitionsForSupplement[["inflammatory_marker"]] = "Inflammatory marker"
definitionsForSupplement[["factor_group"]]	= "Factor group"
definitionsForSupplement[["FVE_Q1"]] = "FVE Q1"
definitionsForSupplement[["FVE_median"]] = "FVE median"
definitionsForSupplement[["FVE_Q3"]] = "FVE Q3"
definitionsForSupplement[["chrom"]]   = "Chromosome"
definitionsForSupplement[["CHROM"]]   = "Chromosome"
definitionsForSupplement[["chr_haplotype"]]   = "Chromosome"
definitionsForSupplement[["start_haplotype"]] = "Haplotype block start"
definitionsForSupplement[["end_haplotype"]] = "Haplotype block end"
definitionsForSupplement[['snp']] = 'SNP' 
definitionsForSupplement[['chr']] = 'Chromosome'
definitionsForSupplement[['pos']] = 'Position'
definitionsForSupplement[['ref']] = 'Reference allele'
definitionsForSupplement[['alt']] = 'Alternative allele'
definitionsForSupplement[['region']] = 'Region' 
definitionsForSupplement[['factor']] = 'Factor' 
definitionsForSupplement[['pval']] = "P-value" 
definitionsForSupplement[['tval']] = 'T-value' 
definitionsForSupplement[['beta']] = 'Beta' 
definitionsForSupplement[['statistic']] = 'T-value' 
definitionsForSupplement[['coefficient']] = 'Coefficient' 
definitionsForSupplement[['coef']] = 'Coefficient' 
definitionsForSupplement[['padj']] = 'padj' 
definitionsForSupplement[['gene_set']] = "Gene set"
definitionsForSupplement[['library']] = "Library"
definitionsForSupplement[['cytoData.good.log2']] = 'cytokine'
definitionsForSupplement[['circMarkerData.log2']] = 'inflaMark'
definitionsForSupplement[['cytoData.good.fc']] = 'cytokine'
definitionsForSupplement[['circMarkerData.fc']] = 'inflaMark'
definitionsForSupplement[['cytoData.good.fcCorr']] = 'cytokine'
definitionsForSupplement[['circMarkerData.fcCorr']] = 'inflaMark'
definitionsForSupplement[['cytoData.good.log2_LONG']] = 'Cytokine'
definitionsForSupplement[['circMarkerData.log2_LONG']] = 'Inflammatory marker'
definitionsForSupplement[['cytoData.good.fc_LONG']] = 'Cytokine'
definitionsForSupplement[['circMarkerData.fc_LONG']] = 'Inflammatory marker'
definitionsForSupplement[['cytoData.good.fcCorr_LONG']] = 'Cytokine'
definitionsForSupplement[['circMarkerData.fcCorr_LONG']] = 'Inflammatory marker'
definitionsForSupplement[['V1']] = 'd0'
definitionsForSupplement[['V2']] = 'd14FC'
definitionsForSupplement[['V3']] = 'd90FC'
definitionsForSupplement[['ATACseq']] = 'ATACseq'
definitionsForSupplement[['SNP']] = 'SNP'
definitionsForSupplement[['ATAC']] = 'ATACseq'

dataTypeNameMapSupplement = list()
dataTypeNameMapSupplement[['cytoData.good.log2']] = 'cytokine'
dataTypeNameMapSupplement[['circMarkerData.log2']] = 'inflaMark'
dataTypeNameMapSupplement[['cytoData.good.fc']] = 'cytokine'
dataTypeNameMapSupplement[['circMarkerData.fc']] = 'inflaMark'
dataTypeNameMapSupplement[['cytoData.good.fcCorr']] = 'cytokine'
dataTypeNameMapSupplement[['circMarkerData.fcCorr']] = 'inflaMark'

tableAdditionToFileNames = list()
tableAdditionToFileNames[["TableS5"]] = "GWAS_top_hits"
tableAdditionToFileNames[["TableS6"]] = "GWAS_pathways"
tableAdditionToFileNames[["TableS7"]] = "explained_variance_and_overlap"

orderToPlotCytokines = readLines('orderToPlotCytokines_newFormat.txt', encoding="UTF-8")
orderToPlotInflaMark = readLines('orderToPlotInflaMarkers_newFormat.txt', encoding="UTF-8")


changePlotRasterization <- function(x,rasterizePlot=T){
  library(ggrastr)
  if (rasterizePlot){
    return(rasterise(x, dpi = 150))
  } else {
    return(x)
  }
}


loadAllAtacSeqAssociationResults <- function(selFolder){
  #p.value.V1_CORR_CYTO.C.albicans.yeast_24h_PBMC_IL.10_good
  #p.value.CYTO.S.aureus_7d_PBMC_IFNg_good
  #Check all the atacSeq files that meet the requirements
  patternSel = paste0('de\\_results\\_.*(\\_|\\.)',timepointSel,'(\\_|\\.).*\\.csv\\.gz$')

  allAtacResFiles = list.files(selFolder, include.dirs = F, pattern = patternSel, recursive = F)
  
  #Iterate over all files and load the relevant data, combine it into a total df
  combinedResults = NULL
  for (fileSel in allAtacResFiles){
    print(fileSel)
    variableSel = gsub('.*\\.visit\\_time.*(CYTO|CM)\\.',
                       '',
                       fileSel)
    
    variableSel = gsub('\\.csv\\.gz$',
                       '',
                       variableSel)
    
    ##Skip if the variables are in the list of cytokines to exclude
    if (grepl(pattern = paste0('(',paste(cytosToExclude,collapse='|'),')'),variableSel)){
      next
    }
    
    print('reading results')
    partialResults = fread(file = file.path(selFolder,fileSel), sep = ',', data.table = F)
    
    #Quick check for NA or INF
    pValColSel = colnames(partialResults)[grep(colnames(partialResults), pattern = paste0('p\\.value\\..*\\.',variableSel))]
    stopifnot(abs(max(partialResults[,pValColSel]))!=abs(min(partialResults[,pValColSel])))
    
    print('merging results')
    colnames( partialResults) = gsub(paste0('^(p.value|t|Coef|Res)\\..*'),'\\1',x = colnames( partialResults))
    partialResults$variable = variableSel
    
    #Now add characterization
    mergedResults = dplyr::left_join(partialResults,peakAnnotations, by = c("genes" = "peak_id"))
    
    if (length(partialResults)>0){
      if (is.null(combinedResults)){
        combinedResults = mergedResults
      }else{
        combinedResults = rbind.data.frame(combinedResults, mergedResults)
      }
    }else{
      break()
    }
    
  }
  
  ###Filter out specific markers/cytokines depending on the specific datatype
  if (dataTypeSel=="cytoData.good.fc"){
    combinedResults = combinedResults[which(gsub('_good','',combinedResults$variable) %in% 
                                              gsub('.fc','',cytokineGroups$cytoData.good.fc$innate_nonspecific_24h_wo_LAC)),]
  }
  return(combinedResults)
}

mapPeakNamesToLocation <- function(peakNames,filePeakAnnotations){
  library(data.table)
  #Load peak annotations
  peakAnnotations = data.table::fread(filePeakAnnotations, data.table = F)
  mappingList = as.list(paste(peakAnnotations[,'chr'], peakAnnotations[,'start'], peakAnnotations[,'end'], sep='_'))
  names(mappingList) = peakAnnotations[,"peak_id"]
  return(unlist(mappingList[peakNames]))
}
  

fixCellCountNames <- function(allMetaDataLukas.input){
  #==Cell Count data percentages PBMC==
  colIndices.PBMC.perc  = grep(pattern = '^PBMC_PERC.*', x = colnames(allMetaDataLukas.input))
  if (length(colIndices.PBMC.perc)>0){
    cellCountData.PBMC.perc = allMetaDataLukas.input[,colIndices.PBMC.perc]
    colnames(cellCountData.PBMC.perc) = gsub('.*:','', x = colnames(cellCountData.PBMC.perc))
    colnames(cellCountData.PBMC.perc) = gsub('/','.', x = colnames(cellCountData.PBMC.perc))
    colnames(cellCountData.PBMC.perc) = paste(colnames(cellCountData.PBMC.perc),'pbmc.perc',sep='.')
  } else {
    cellCountData.PBMC.perc=NULL
  }
  
  #==Cell Count data percentages WB==
  colIndices.WB.perc = grep(pattern = '^WB_PERC.*', x = colnames(allMetaDataLukas.input))
  if (length(colIndices.WB.perc)>0){
    cellCountData.WB.perc = allMetaDataLukas.input[,colIndices.WB.perc]
    colnames(cellCountData.WB.perc) = gsub('.*:','', x = colnames(cellCountData.WB.perc))
    colnames(cellCountData.WB.perc) = gsub('/','.', x = colnames(cellCountData.WB.perc))
    colnames(cellCountData.WB.perc) = paste(colnames(cellCountData.WB.perc),'wb.perc',sep='.')
  } else {
    cellCountData.WB.perc=NULL
  }
  
  #==Cell Count data absolute counts ===
  colIndices.WB.abs =  grep(pattern = '^WB_PER_ML.*', x = colnames(allMetaDataLukas.input))
  if (length(colIndices.WB.abs)>0){
    cellCountData.WB.abs = allMetaDataLukas.input[,colIndices.WB.abs]
    colnames(cellCountData.WB.abs) = gsub('.*:','', x = colnames(cellCountData.WB.abs))
    colnames(cellCountData.WB.abs) = gsub('/','.', x = colnames(cellCountData.WB.abs))
    colnames(cellCountData.WB.abs) = paste(colnames(cellCountData.WB.abs),'wb.abs',sep='.')
  }else{
    cellCountData.WB.abs=NULL
  }
  
  return(list('cellCountData.PBMC.perc'=cellCountData.PBMC.perc,
              'cellCountData.WB.perc'=cellCountData.WB.perc,
              'cellCountData.WB.abs'=cellCountData.WB.abs))
}


loadLukasMetaDataPbmc <- function(allMetaDataPath){
  #Load the data from Lukas' file with all the cytokines and other metadata
  library(data.table)
  library(Rfast)
  allMetaDataLukas = fread(allMetaDataPath, data.table = F)
  rownames(allMetaDataLukas) = allMetaDataLukas$"SAMPLE:ID"
  
  allMetaDataLukas.pbmc = allMetaDataLukas[grep(rownames(allMetaDataLukas),pattern = '\\_PBMC'),]
  idsPresentMetaData = gsub('^300','',rownames(allMetaDataLukas.pbmc))
  idsPresentMetaData = gsub('\\_PBMC.*','',idsPresentMetaData)
  rownames(allMetaDataLukas.pbmc) = idsPresentMetaData
  #=====CYTOKINES=====

  extractDataSubsetFromDataframe <- function(dataFrameSel, patternSel){
    dataFrameSel = dataFrameSel[,grep(pattern = patternSel, x = colnames(dataFrameSel))]
    colnames(dataFrameSel) = gsub('.*:','', x = colnames(dataFrameSel))
    colnames(dataFrameSel) = gsub('(\\_good|\\_ok)','', x = colnames(dataFrameSel))
    return(dataFrameSel)
  }
  
  #==Cytokine data==
  cytoData.good.log2 = extractDataSubsetFromDataframe(allMetaDataLukas.pbmc, '^CYTO.*\\_good$')
  #Remove some cytokines that need to be excluded
  cytoData.good.log2 = cytoData.good.log2[,which(!(colnames(cytoData.good.log2) %in% cytosToExclude))]

  #==Cytokine fold change==
  cytoData.good.fc = extractDataSubsetFromDataframe(allMetaDataLukas.pbmc, '^LFC\\_CYTO.*\\_good$')
  colnames(cytoData.good.fc) = paste0(colnames(cytoData.good.fc),'.fc')
  #Remove some cytokines that need to be excluded
  colsToKeepFc = grep(pattern = paste0('(',paste(cytosToExclude,collapse='|'),').*'),colnames(cytoData.good.fc), invert = T)
  cytoData.good.fc = cytoData.good.fc[,colsToKeepFc]

  #Also load the fold-change data with the cell-counts changes and the changes in visit time
  #regressed out beforehand
  cytoData.good.fcCorr =  extractDataSubsetFromDataframe(allMetaDataLukas.pbmc, '^LFC\\_CORR\\_CYTO.*\\_good$')
  colnames(cytoData.good.fcCorr) = paste0(colnames(cytoData.good.fcCorr),'.fcCorr')
  #Remove some cytokines that need to be excluded
  colsToKeepFc = grep(pattern = paste0('(',paste(cytosToExclude,collapse='|'),').*'),colnames(cytoData.good.fcCorr), invert = T)
  cytoData.good.fcCorr = cytoData.good.fcCorr[,colsToKeepFc]
  
  
  
  #=====CIRCULATING MARKERS OF INFLAMMATION======
  #==Circulating markers of inflammation==
  circMarkerData.log2 = allMetaDataLukas.pbmc[,grep(pattern = '^CM.*', x = colnames(allMetaDataLukas.pbmc))]
  colnames(circMarkerData.log2) = gsub('.*:','', x = colnames(circMarkerData.log2))
  colnames(circMarkerData.log2) = gsub('/','.', x = colnames(circMarkerData.log2))
  
  #==Circulating markers fold-changes
  circMarkerData.fc = allMetaDataLukas.pbmc[,grep(pattern = '^LFC_CM.*', x = colnames(allMetaDataLukas.pbmc))]
  colnames(circMarkerData.fc) = gsub('.*:','', x = colnames(circMarkerData.fc))
  colnames(circMarkerData.fc) = gsub('/','.', x = colnames(circMarkerData.fc))
  colnames(circMarkerData.fc) = paste(colnames(circMarkerData.fc),'fc',sep='.')
  
  #Also load the fold-change data with the cell-counts changes and the changes in visit time
  #regressed out beforehand
  circMarkerData.fcCorr = allMetaDataLukas.pbmc[,grep(pattern = '^LFC_CORR_CM.*', x = colnames(allMetaDataLukas.pbmc))]
  colnames(circMarkerData.fcCorr) = gsub('.*:','', x = colnames(circMarkerData.fcCorr))
  colnames(circMarkerData.fcCorr) = gsub('/','.', x = colnames(circMarkerData.fcCorr))
  colnames(circMarkerData.fcCorr) = paste(colnames(circMarkerData.fcCorr),'fcCorr',sep='.')
  
  
  #Calculate fold changes of the circulating mediators as a check
  circMarkerData.fc.check = circMarkerData.log2
  circMarkerData.fc.check[] = NA
  refTp = 'V1'
  refRows = grep(paste0('.*',refTp,'$'),x = rownames(circMarkerData.log2))
  for (otherTp in c('V2','V3')){
    selRows = grep(paste0('.*',otherTp,'$'),x = rownames(circMarkerData.log2))
    circMarkerData.fc.check[selRows,] = circMarkerData.log2[selRows,]-circMarkerData.log2[refRows,]
    
    stopifnot(all.equal(gsub(x = rownames(circMarkerData.log2)[refRows], pattern = '\\_V[1-3]$', replacement = ''),
                        gsub(x = rownames(circMarkerData.log2)[selRows], pattern = '\\_V[1-3]$', replacement = '')))
  }
  colnames(circMarkerData.fc.check) = paste(colnames(circMarkerData.fc.check),'fc',sep='.')
  
  stopifnot(all.equal(circMarkerData.fc.check,circMarkerData.fc))
  
  #==========CELL COUNTS=======
  cellCountDataFixedList = fixCellCountNames(allMetaDataLukas.pbmc)
  cellCountData.PBMC.perc = cellCountDataFixedList[['cellCountData.PBMC.perc']]
  cellCountData.WB.perc =cellCountDataFixedList[['cellCountData.WB.perc']]
  cellCountData.WB.abs =cellCountDataFixedList[['cellCountData.WB.abs']]
  
  #===Cell count/percentages fold changes with V1===
  cellCountData.PBMC.perc.fc = cellCountData.PBMC.perc
  cellCountData.PBMC.perc.fc[] = NA
  cellCountData.WB.abs.fc = cellCountData.WB.abs
  cellCountData.WB.abs.fc[] = NA
  for (tpSel in c("V3","V2","V1")){
    tpSelIndices = grep(pattern = paste0('*\\_',tpSel,'$'), x = rownames(cellCountData.PBMC.perc.fc))
    tp1Indices = grep(pattern = paste0('*\\_V1$'), x = rownames(cellCountData.PBMC.perc.fc))
    
    stopifnot(all.equal(gsub(x = rownames(cellCountData.PBMC.perc.fc)[tpSelIndices], pattern = '\\_V[1-3]$', replacement = ''),
                        gsub(x = rownames(cellCountData.PBMC.perc.fc)[tp1Indices], pattern = '\\_V[1-3]$', replacement = '')))
    stopifnot(all.equal(gsub(x = rownames(cellCountData.WB.abs.fc)[tpSelIndices], pattern = '\\_V[1-3]$', replacement = ''),
                        gsub(x = rownames(cellCountData.WB.abs.fc)[tp1Indices], pattern = '\\_V[1-3]$', replacement = '')))
    
    if (tpSel %in%c('V3','V2')){
      cellCountData.PBMC.perc.fc[tpSelIndices,] = log2(cellCountData.PBMC.perc[tpSelIndices,]/cellCountData.PBMC.perc[tp1Indices,])
      cellCountData.WB.abs.fc[tpSelIndices,] = log2(cellCountData.WB.abs[tpSelIndices,]/cellCountData.WB.abs[tp1Indices,])
    }else{
      cellCountData.PBMC.perc.fc[tpSelIndices,] = NA
      cellCountData.WB.abs.fc[tpSelIndices,] = NA
    }
  }
  colnames(cellCountData.PBMC.perc.fc) = gsub('\\.perc$','.fc',colnames(cellCountData.PBMC.perc.fc))
  colnames(cellCountData.WB.abs.fc) = gsub('\\.abs$','.fc',colnames(cellCountData.WB.abs.fc))
  
  
  
  
  
  #=========TIME OF VISIT=========
  #===Time of visit diff with timepoint 1====
  allMetaDataLukas.pbmc$`SAMPLE:VISIT_TIME_REAL.diff`=allMetaDataLukas.pbmc$`SAMPLE:VISIT_TIME_REAL`
  for (tpSel in c('V3','V2',"V1")){
    tpSelIndices = grep(pattern = paste0('*\\_',tpSel,'$'), x = rownames(allMetaDataLukas.pbmc))
    tp1Indices = grep(pattern = paste0('*\\_V1$'), x = rownames(allMetaDataLukas.pbmc))
    stopifnot(all.equal(gsub(x = rownames(allMetaDataLukas.pbmc)[tpSelIndices], pattern = '\\_V[1-3]$', replacement = ''),
                        gsub(x = rownames(allMetaDataLukas.pbmc)[tp1Indices], pattern = '\\_V[1-3]$', replacement = '')))
    if (tpSel %in%c('V3','V2')){
      allMetaDataLukas.pbmc[tpSelIndices,"SAMPLE:VISIT_TIME_REAL.diff"] = allMetaDataLukas.pbmc[tpSelIndices,"SAMPLE:VISIT_TIME_REAL.diff"]-allMetaDataLukas.pbmc[tp1Indices,"SAMPLE:VISIT_TIME_REAL.diff"]
    }else{
      allMetaDataLukas.pbmc[tpSelIndices,"SAMPLE:VISIT_TIME_REAL.diff"] = NA
    }
  }
  
  
  
  
  #=======SELECT OTHER PARAMETERS METADATA========
  metaData = allMetaDataLukas.pbmc[,c("RUN:TSS_ENRICHMENT","LAB:BATCH","DONOR:AGE",
                                      "DONOR:SEX","DONOR:BMI","SAMPLE:VISIT","SAMPLE:VISIT_DATE_2PI_COS", 
                                      "SAMPLE:VISIT_DATE_2PI_SIN", "SAMPLE:VISIT_DAYS_FROM_20170101",
                                      "SAMPLE:VISIT_TIME_REAL","SAMPLE:VISIT_TIME","DONOR:IC_TIME","DONOR:IC_TIME_REAL",
                                      "DONOR:IC_EVENING","SAMPLE:VISIT_TIME_REAL.diff",
                                      "SAMPLE:EXCLUSION",
                                      "SAMPLE:alcoholInLast24h",
                                      "DONOR:SNP_OUTLIER","DONOR:oralContraceptivesIncludingMen",
                                      "thm.adaptive_MTB_7d_V3_FC1.2_responder",
                                      "thm.heterologous_nonspecific_7d_V3_FC1.2_responder",
                                      "thm.innate_nonspecific_24h_wo_LAC_IL10_IL1ra_V3_FC1.2_responder")]
  
  colnames(metaData) = gsub('.*:','', x = colnames(metaData))
  colnames(metaData)[which(colnames(metaData)=="VISIT")] = "timePoint"
  metaData$person = gsub('\\_V[1-3]$','',rownames(metaData))
  metaData[which( !(metaData$SEX %in% c('M','F'))),'SEX'] = NA
  
  metaData$oralContraceptivesIncludingMen = as.numeric(metaData$oralContraceptivesIncludingMen)
  
  dataToReturn = list("circMarkerData.log2"=circMarkerData.log2,"circMarkerData.fc"=circMarkerData.fc,"circMarkerData.fcCorr"=circMarkerData.fcCorr,
                      "cytoData.good.log2"=cytoData.good.log2,'cytoData.good.fc'=cytoData.good.fc,"cytoData.good.fcCorr"=cytoData.good.fcCorr,
                      "cellCountData.PBMC.perc"=cellCountData.PBMC.perc,"cellCountData.WB.perc"=cellCountData.WB.perc,
                      "cellCountData.WB.abs"=cellCountData.WB.abs,
                      "cellCountData.PBMC.perc.fc" = cellCountData.PBMC.perc.fc, "cellCountData.WB.abs.fc" = cellCountData.WB.abs.fc,
                      "metaData"=metaData)
  #Check that these df have no overlapping column names
  allColnames = lapply(dataToReturn, function(x) colnames(x))
  stopifnot(max(table(unlist(allColnames)))==1)
  allRelevantData = base::cbind.data.frame(dataToReturn)
  #Remove wrong names
  allNamesToReplace = names(dataToReturn)
  allNamesToReplace = allNamesToReplace[order(nchar(allNamesToReplace), decreasing = T)]
  for (nameSelToRep in allNamesToReplace){
    colnames(allRelevantData) = gsub(pattern = paste0(nameSelToRep,'.'),replacement = '',colnames(allRelevantData),fixed=T)
  }
  dataToReturn[['allRelevantData']]=allRelevantData
  
  return(dataToReturn)
}





loadCorrectionFactorList <- function(){
  basicCellTypesToCorrectCyto=c('MONO', 'MONO/INTERMEDIATE', 'MONO/NON_CLASSICAL',
                                'T/CD8', 'T/CD4', 'T/CD4/TREG',
                                'B', 'NK', 'NKT',
                                'BASO', 'NEUTRO')
  basicCellTypesToCorrectCyto = gsub('/','.',basicCellTypesToCorrectCyto, fixed=T)
  #For stimulated cytokines we use cell type percentages in PBMC
  correctionCytokines.abs = c("SEX",
                              "AGE",
                              "VISIT_TIME_REAL",
                              paste0(basicCellTypesToCorrectCyto,'.pbmc.perc'))
  
  #For the fold changes of cytos we use the fold changes of the cell types in PBMC and the differences in visit time
  correctionCytokines.fc = c("SEX",
                             "AGE",
                             "VISIT_TIME_REAL.diff",
                             paste0(basicCellTypesToCorrectCyto,'.pbmc.fc'))
  
  correctionCytokines.fcCorr = c("SEX",
                             "AGE")
  
  #For the absolute circulating mediators we use the absolute concentrations in Whole Blood
  basicCellTypesToCorrectCircMed=  c( 'MONO', 'MONO/INTERMEDIATE', 'MONO/NON_CLASSICAL',
                                      'T/CD8', 'T/CD4', 'T/CD4/TREG',
                                      'B', 'NK', 'NKT',
                                      'NEUTRO')
  basicCellTypesToCorrectCircMed = gsub('/','.',basicCellTypesToCorrectCircMed, fixed=T)
  correctionCircMed.abs = c("SEX",
                            "AGE",
                            "BMI",
                            "oralContraceptivesIncludingMen",
                            "VISIT_TIME_REAL",
                            paste0(basicCellTypesToCorrectCircMed,'.wb.abs'))
  
  #For the fold-changes in circ. mediators we use the fold-changes of cell counts in whole blood
  correctionCircMed.fc = c("SEX",
                           "AGE",
                           "BMI",
                           "oralContraceptivesIncludingMen",
                           "VISIT_TIME_REAL.diff",
                           paste0(basicCellTypesToCorrectCircMed,'.wb.fc'))
  
  correctionCircMed.fcCorr = c("SEX",
                           "AGE",
                           "BMI",
                           "oralContraceptivesIncludingMen")
  
  correctionFactorList = list("cytoData.good.log2"=correctionCytokines.abs,
                              "cytoData.good.fc"=correctionCytokines.fc,
                              "cytoData.good.fcCorr"= correctionCytokines.fcCorr,
                              "circMarkerData.log2"=correctionCircMed.abs,
                              "circMarkerData.fc"=correctionCircMed.fc,
                              "circMarkerData.fcCorr"=correctionCircMed.fcCorr)
  
  
  
}




plotTopPathwaysPvals <- function(selPvals.df, colSel, titleSel=NULL){
  #selPvals.df - df with one column that contains the pvalues of interest
  #colSel - selected column name with the pvalues of interest
  selColumnPvals = selPvals.df[,colSel,drop=F]
  selColumnPvals[[paste0(colSel,'.adjusted')]] = p.adjust(selColumnPvals[,1],method = 'fdr')
  selColumnPvals = selColumnPvals[order(selColumnPvals[,1], decreasing = F),]
  selColumnPvals[,'minLog10'] = -log10(selColumnPvals[,colSel])
  selColumnPvals$pathways = rownames(selColumnPvals)
  
  #Select the top
  topToPlot = 20
  if (nrow(selColumnPvals)<20){
    topToPlot =nrow(selColumnPvals)
  }
  selColumnPvalsToPlot = selColumnPvals[c(1:topToPlot),]
  selColumnPvalsToPlot$pathways = factor(selColumnPvalsToPlot$pathways, levels = rev(selColumnPvalsToPlot$pathways))
  #determine the position of the FDR cutoff
  minAdjPval = min(selColumnPvalsToPlot[,paste0(colSel,'.adjusted')],na.rm=T)
  if (minAdjPval>.05){
    vLinePos = 1.05* max(selColumnPvalsToPlot[,'minLog10'])
  }else{
    #find first significant pval
    smallerThanCutoff = max(which(selColumnPvals[,paste0(colSel,'.adjusted')]<=.05))
    #biggerThanCutoff =  which(selColumnPvalsToPlot[,paste0(colSel,'.adjusted')]>.05)
    vLinePos = mean(c(selColumnPvals[[smallerThanCutoff,'minLog10']],selColumnPvals[[smallerThanCutoff+1,'minLog10']]))
  }
  #Now make the plot
  p1 = ggplot(selColumnPvalsToPlot, aes(x = pathways, y = minLog10)) +
    geom_bar(stat="identity", position="identity", fill = '#C54D52') +
    theme_bw() + theme(

      axis.text = element_text(size = 12),
      text = element_text(family = "Helvetica"),
      legend.position = "none",
      axis.title.x = element_text(size=12),
      axis.title.y=element_blank(),
      panel.grid.minor.y = element_blank(),
      panel.grid.major.y = element_blank()
      
      
    ) + coord_flip() + 
    geom_hline(yintercept = vLinePos, linetype = "dashed") + ylab('-log10(p-value)')
  
  if (!is.null(titleSel)){
    p1 = p1 + ggtitle(titleSel) + theme(plot.title = element_text(size = 5))
  }
  
  return(p1)
}





processPathwayNames <- function(origNames, addGo=F, maxLengthChar=38, maxNrow = 2){
  ##Process pathways names to make them fit better
  #Split them over two lines if I want with "maxNrow" and give a maximum length with "maxLengthChar".
  #The maximum nr of characters in "maxLengthChar" is basically multiplied by maxNrow to give the actual
  #maximum nr of characters that will be displayed.
  #origNames - the original pathway names to be altered
  #addGo - add "(GO)" at the end of a pathway if it comes from GO
  #maxLengthChar - maximum nr of characters on a single line
  #maxNrow - the number of lines used per pathway if they exceed the "maxLengthChar"
  
  #remove the GO part
  if (addGo){
    processedNames = gsub(' \\(GO\\:[0-9]{5,7}\\)$','(GO)',origNames)
  }else{
    processedNames = gsub(' \\(GO\\:[0-9]{5,7}\\)$','',origNames) 
  }
  #Replace some common parts of GO names
  processedNames = gsub('negative regulation of ','neg. ',processedNames)
  processedNames = gsub('positive regulation of ','pos. ',processedNames)
  processedNames = gsub('regulation of ','reg. ',processedNames)
  processedNames = gsub('cellular response to ','cell resp. ',processedNames)
  processedNames = gsub('activation of ','',processedNames)
  processedNames = gsub('response to ','resp. ',processedNames)
  processedNames = gsub(' pathway$','',processedNames)
  processedNames = gsub(' process$',' pr.',processedNames)
  
  processedNames = gsub("Fc gamma R","FcγR",processedNames)
  processedNames = gsub("tumor necrosis factor","TNF",processedNames)
  processedNames = gsub("epidermal growth factor","EGF",processedNames)
  processedNames = gsub("Staphylococcus aureus","S. aureus",processedNames)
  processedNames = gsub("Herpes simplex virus 1 infection","HSV-1",processedNames)
  processedNames = gsub("Kaposi sarcoma-associated herpesvirus","KSHV",processedNames)
  processedNames = gsub("polymerase II","pol2",processedNames)
  processedNames = gsub("fibroblast growth factor","FGF",processedNames)
  
  strReverse <- function(x){
    sapply(lapply(strsplit(x, NULL), rev), paste, collapse="")
    }
  
  #Now start going into the maximum number of characters
  allProcessedNamesCheckDouble = c()
  for (indexSel in c(1:length(processedNames))){
    oldName = processedNames[[indexSel]]
    
    #Remove anything between brackets
    oldName = gsub("\\s*\\([^\\)]+\\)","",oldName)
    
    if (nchar(processedNames[[indexSel]])>maxLengthChar){
      #First split into the number of lines allowed if the number of characters exceeds the maximum
      #Kind of hacky way to use a regex. My regex skills were not up to the task, and stackoverflow only had the reverse answer (starting from the end)
      #So, yeah, I hacked it by reversing the string
      #https://stackoverflow.com/questions/53337711/how-to-insert-a-character-every-n-characters-from-end-of-string
      #newName = strReverse(gsub(paste0('(?s)(?=(?:.{',maxLengthChar,'})+$)'), strReverse("\n"), strReverse(oldName), perl = TRUE))
      #The function below only works for two lines max so far
      newName = gsub(paste0('(.{1,',maxLengthChar,'})\\s(.*)'),'\\1->\n\\2',oldName, perl = TRUE)
    
      #After that remove the last part of the name if it is too long
      #Multiply the number to get the real max, and add 2 for each /n
      realMaxNr = maxLengthChar * maxNrow + (2*(maxNrow-2))
      if (endsWith(oldName, suffix = '(GO)')){
        newName = substr(newName,start = 1, stop=realMaxNr-2)
        newName = paste0(newName,'(GO)')
      }else{
        newName = substr(newName,start = 1, stop=realMaxNr)
      }
      
    }else{
      newName = oldName
    }
    
    while (newName %in% allProcessedNamesCheckDouble){
      newName = paste0(newName,'.2')
    }
    
    processedNames[[indexSel]] = newName
    allProcessedNamesCheckDouble = c(allProcessedNamesCheckDouble,newName)
  }
  return(processedNames)
}




loadCountsPbmcLukas <- function(fileWithAtacData){
  counts = fread(input = fileWithAtacData, data.table = F)
  rownames(counts) = counts[,1]
  counts[,1] <- NULL
  #Check for which indiv we have counts
  idsPresent = gsub('^300','',colnames(counts))
  idsPresent = gsub('\\_PBMC.*','',idsPresent)
  stopifnot(max(table(idsPresent))==1)
  oldColNamesAtac = colnames(counts)
  #Rename the columns
  colnames(counts) = idsPresent
  return(counts)
}

loadPathwayDatabasesLukas <- function(folderPathways, patternEnd = '_entrez.gmt'){
  allPathwayFiles = list.files(folderPathways, pattern = paste0(patternEnd,'$'))
  pathwayList = list()
  for (fileSel in allPathwayFiles){
    nameToSave = gsub( patternEnd,'',fileSel)
    pathwayList[[nameToSave]] = list()
    pathwaySetRaw = readLines(file.path(folderPathways,fileSel))
    for (rowSel in pathwaySetRaw){
      splitPathwayInfo = str_split(rowSel, pattern = '\t')[[1]]
      pathwayName = splitPathwayInfo[[1]]
      genes = splitPathwayInfo[c(-1,-2)]

      pathwayList[[nameToSave]][[pathwayName]] = genes
    }
  }
  return(pathwayList)
}


plotTopPathwaysCounts <- function(percentageResults.ordered, pathwaySelName, colorOutline = '#000000'){
  #percentageResults.ordered - should have 3  columns. The first column
  #   should be called "pathway" and has the pathway names. The second
  #   column 'direction' should indicate whether the pathway analysis was ran for
  #   genes up or downregulated ['Up' and 'Down'] (for genetics this 
  #   should always be 'Up'). The last column is called 'percSignif' and 
  #   contains the actual percentages to be shown
  #   The order of the up and downregulated pathways should be the same.
  #   The order in which they are provided in the order in which they
  #   are plotted.
  #pathwaySelName - name of the pathwayset selected (e.g. GO or KEGG)
  
  #Define the limits: if there are just up or just down pathways, than range from 0 to the max/min
  #otherwise go from min to max

  if (length(setdiff(c('Up','Down'),percentageResults.ordered$direction))==0){
    startLim = - max(abs(percentageResults.ordered$percSignif))
    endLim = max(abs(percentageResults.ordered$percSignif))
  }else if (unique(percentageResults.ordered$direction)=='Up'){
    startLim = 0
    endLim = max(abs(percentageResults.ordered$percSignif))
  }else if (unique(percentageResults.ordered$direction)=='Up'){
    startLim = - max(abs(percentageResults.ordered$percSignif))
    endLim = 0
  }

  #Capitalize first letter
  origLevels = levels(percentageResults.ordered$pathway)
  levelsCapped = capitalize(origLevels)
  percentageResults.ordered$pathway = capitalize(percentageResults.ordered$pathway)
  percentageResults.ordered$pathway = factor(percentageResults.ordered$pathway, levels = levelsCapped)
  
  #Now make the plot
  tickSizes = 9
  lineheight=.65
  p = ggplot(percentageResults.ordered, aes(x = pathway, y = percSignif, fill = direction)) +
    geom_bar(stat="identity", position="identity") +
    ylim(startLim, endLim) +
    theme_bw() +
    theme(axis.text.x = element_text(size = tickSizes, lineheight = lineheight, colour = 'black'),
          axis.text.y = element_text(size = tickSizes, lineheight = lineheight, colour = 'black'),
          text = element_text(family = "Helvetica", colour = 'black'),
          legend.position = "none",
          axis.title.x = element_text(size=12, colour = 'black'),
          axis.title.y = element_blank(),
          panel.grid.minor.y = element_blank(),
          panel.grid.major.y = element_blank()
          ) + coord_flip() +  
    scale_fill_manual(values = c('#c44e52','#4c72b0')) +
    ylab('Percentage of markers with P ≤ 0.05')
  
  breaksSel = ggplot_build(p)$layout$panel_params[[1]]$x$breaks
  p = p + scale_y_continuous(breaks = breaksSel,
                             labels = paste0(as.character(abs(breaksSel)), '%'))
  
  p = p + ggtitle(pathwaySelName) + theme(plot.title = element_text(hjust = 0.5, size=12))
  p
  
  return(p)
}



defineSubgroupsCytokines <- function(){
  allRelevantDataList = loadLukasMetaDataPbmc(allMetaDataPath)
  list2env(allRelevantDataList, .GlobalEnv)
  
  ####Define the cytokine groups####
  #--baseline--
  cytokineGroups = list()
  allCytosNames.log2 = colnames(allRelevantDataList$cytoData.good.log2)
  cytokineGroups[['cytoData.good.log2']] = list()
  cytokineGroups[['cytoData.good.log2']][['adaptive_7d']] = allCytosNames.log2[grep(allCytosNames.log2,pattern = '\\_7d\\_')]
  cytokineGroups[['cytoData.good.log2']][['innate_24h_pro_inflammatory']] = allCytosNames.log2[grep("(?=.*24h)(?!.*(IFNg|IL\\.17|IL\\.1ra|IL\\.10|lactate))", allCytosNames.log2, perl=TRUE)]
  cytokineGroups[['cytoData.good.log2']][['innate_24h_anti_inflammatory']] = allCytosNames.log2[grep("(?=.*24h)(?=.*(IL\\.1ra|IL\\.10))", allCytosNames.log2, perl=TRUE)]
  cytokineGroups[['cytoData.good.log2']][['innate_24h_cyto']] = allCytosNames.log2[grep("(?=.*24h)(?!.*(lactate|IFNg|IL\\.17))", allCytosNames.log2, perl=TRUE)]
  cytokineGroups[['cytoData.good.log2']][['all_cyto_noLact']] = c(cytokineGroups[['cytoData.good.log2']][['innate_24h_cyto']],cytokineGroups[['cytoData.good.log2']][['adaptive_7d']])
  #--foldChange--
  allCytosNames.fc = colnames(allRelevantDataList$cytoData.good.fc)
  cytokineGroups[['cytoData.good.fc']] = list()
  cytokineGroups[['cytoData.good.fc']][['adaptive_MTB_7d']] = allCytosNames.fc[grep(allCytosNames.fc,pattern = '^MTB.*\\_7d\\_.*')]
  cytokineGroups[['cytoData.good.fc']][['heterologous_nonspecific_7d']] = allCytosNames.fc[grep("^(?!MTB).*\\_7d.*$", allCytosNames.fc, perl = TRUE)]
  cytokineGroups[['cytoData.good.fc']][['innate_nonspecific_24h_wo_LAC']] = allCytosNames.fc[grep("^(?!MTB).*\\_24h.*(?!(lactate)).*$", allCytosNames.fc, perl = TRUE)]
  cytokineGroups[['cytoData.good.fc']][['innate_nonspecific_24h_wo_LAC_IL10_IL1ra']] = allCytosNames.fc[grep("^(?!(MTB)).*\\_24h(?!.*(IL\\.1ra|IL\\.10|lactate)).*$", allCytosNames.fc, perl = TRUE)]
  cytokineGroups[['cytoData.good.fc']][['all_cyto_noLact']] = allCytosNames.fc[grep(".*(?!(lactate)).*$", allCytosNames.fc, perl = TRUE)]
  
  ####circulating markers####
  #--baseline--
  cytokineGroups[['circMarkerData.log2']] = list()
  cytokineGroups[['circMarkerData.log2']][['all']] = colnames(allRelevantDataList[['circMarkerData.log2']])
  #--foldChange--
  cytokineGroups[['circMarkerData.fc']] = list()
  cytokineGroups[['circMarkerData.fc']][['all']] = colnames(allRelevantDataList[['circMarkerData.fc']])
  return(cytokineGroups)
}



mapLongCytoToShortCytoName <- function(longNames){
  longToShortMap = c('Ca','LPS','Sa','Mt')
  names(longToShortMap) = c("C.albicans.yeast","LPS.100ng","S.aureus","MTB")
  return(unlist(lapply(longNames, function(x) longToShortMap[[x]] )))
}



createMappingList <- function(loadedData){
  #Mapping list with values first column
  mappingList = list()
  for (rowIndex in 1:nrow(loadedData)){
    valueSel = loadedData[rowIndex,1]
    for (colSel in 2:ncol(loadedData)){
      keySel = loadedData[rowIndex,colSel]
      if (!is.na(keySel)){
        mappingList[[keySel]] = valueSel
      }
    }
  }
  return(mappingList)
}

mapCytokineNamesFile = file.path('StimulusToPrettyNames.xlsx')
sheetName = 'cytokine'
loadedData = openxlsx::read.xlsx(xlsxFile = mapCytokineNamesFile, sheet = sheetName, colNames=F)
mappingListCyto = createMappingList(loadedData)

sheetName = 'circMed'
loadedData = openxlsx::read.xlsx(xlsxFile = mapCytokineNamesFile, sheet = sheetName, colNames=F)
mappingListCircMed = createMappingList(loadedData)





mapCytoNamesForSupplement <- function(origNames,dataTypeSel){
  if (dataTypeSel %in% c('cytoData.good.log2','cytoData.good.fc','cytoData.good.fcCorr','cytoData.good')){
    shortStimSel = mapLongCytoToShortCytoName(str_split_fixed(origNames, pattern = '_', 2)[,1])
    cytoMapped = mappingListCyto[str_split_fixed(origNames, pattern = '_', 4)[,4]]

    cytoTotal = paste0(cytoMapped, ' (', shortStimSel,', ',  str_split_fixed(origNames, pattern = '_', 4)[,2],')')
  }else if (dataTypeSel %in% c('circMarkerData.log2','circMarkerData.fc','circMarkerData.fcCorr','circMarkerData')){
    cytoMapped = mappingListCircMed[str_split_fixed(origNames, pattern = '_', 2)[,2]]
    cytoTotal = unlist(cytoMapped)
  }
  
  return(cytoTotal)
}


loadCombinedAtacAssociationResults <- function(allRelevantDataList, dataTypeSel,timepointSel, folderWithBaselineAtacResults=NULL, folderWithTrainingAtacResults=NULL,
                                               mapFactorNames=T, mapPeakNames=T, covariateAddToGrep=''){
  ##Load the results of the Limma-VOOM atacseq results
  ##Specifically looking at the associations between ATACseq and cytokine/inflammatory marker production and fold-changes
  #allRelevantDataList - list with all dataframes
  #dataTypeSel - type of data selected "cytoData.good.log2","cytoData.good.fc","circMarkerData.log2","circMarkerData.fc"
  #timepointSel - V1, V2, V3
  #folderWithBaselineAtacResults - if we are going to look at levels at d0
  #folderWithTrainingAtacResults - if we want to look at fold-changes

  #Select data
  dataSel = allRelevantDataList[[dataTypeSel]]
  selVariables = colnames(dataSel)
  
  if (startsWith(x = dataTypeSel, prefix = 'cytoData')){
    selData = 'cytokines'
  }else{
    selData = 'CM'
  }
  
  if (grepl(x = dataTypeSel, pattern = '(\\.fc|\\.fcCorr)$')){
    selFolder = file.path(folderWithTrainingAtacResults,selData)
  }else{
    selFolder = file.path(folderWithBaselineAtacResults,selData)
  }
  
  #Check all the atacSeq files that meet the requirements
  patternSel = paste0('.*','(\\_|\\.)',timepointSel,'(\\_|\\.).*',covariateAddToGrep,'.*\\.csv\\.gz$')
  allAtacResFiles = list.files(selFolder, include.dirs = F, pattern = patternSel)
  
  #Iterate over all files and load the relevant data, combine it into a total df
  combinedGenes_pval = NULL
  combinedGenes_t = NULL
  combinedResults = NULL
  fileSel = allAtacResFiles[[1]]
  fullTable.long = NULL
  for (fileSel in allAtacResFiles){
    print(fileSel)
    variableSel = gsub('.*\\.visit\\_time.*(CYTO|CM)\\.',
                       '',
                       fileSel)
    
    variableSel = gsub('\\.csv\\.gz$',
                       '',
                       variableSel)
    
    ##Skip if the variables are in the list of cytokines to exclude
    if (grepl(pattern = paste0('(',paste(cytosToExclude,collapse='|'),')'),variableSel)){
      next
    }
    
    print('reading results')
    partialResults = fread(file = file.path(selFolder,fileSel), sep = ',', data.table = F)
    colnames(partialResults)[which(colnames(partialResults)=='V1')] = definitionsForSupplement[['region']]
    
    if (mapPeakNames){
    ##Map region names to chromPos
    partialResults[[definitionsForSupplement[['region']]]] = mapPeakNamesToLocation(partialResults[[definitionsForSupplement[['region']]]],
                                                                                    filePeakAnnotations = peakInfoFile)
    }
    
    #Quick check for NA or INF
    pValColSel = colnames(partialResults)[grep(colnames(partialResults), pattern = paste0('p\\.value\\..*\\.',variableSel))]
    stopifnot(abs(max(partialResults[,pValColSel]))!=abs(min(partialResults[,pValColSel])))
    tValColSel = colnames(partialResults)[grep(colnames(partialResults), pattern = paste0('t\\..*\\.',variableSel))]
    coefColSel = colnames(partialResults)[grep(colnames(partialResults), pattern = paste0('Coef\\..*\\.',variableSel))]
    
    ##Also combine it into a long format
    partialResults.long = partialResults[,c(definitionsForSupplement[['region']],coefColSel,tValColSel,pValColSel)]
    colnames(partialResults.long)[unlist(lapply(c(coefColSel,tValColSel,pValColSel),  function(x) which(colnames(partialResults.long)==x)))] = 
      c(definitionsForSupplement[['coefficient']],
        definitionsForSupplement[['tval']],
        definitionsForSupplement[['pval']])
    #add variable name
    partialResults.long = tibble::add_column(partialResults.long, !!(definitionsForSupplement[['factor']]) := variableSel, .after = definitionsForSupplement[['region']])
    
    if (is.null(partialResults.long)){
      fullTable.long = partialResults.long
    }else{
      fullTable.long = rbind.data.frame(fullTable.long, partialResults.long)
    }
  }

  if (mapFactorNames){
  ##Map Cyto names
  cytoPart = gsub('\\_good$','',unlist(fullTable.long[definitionsForSupplement[['factor']]]))
  cytoPart.map = mapCytoNamesForSupplement(cytoPart,dataTypeSel)
  fullTable.long[definitionsForSupplement[['factor']]] = cytoPart.map
  }else{
    fullTable.long[definitionsForSupplement[['factor']]] = gsub('\\_good$','',unlist(fullTable.long[definitionsForSupplement[['factor']]]))
  }
  return(fullTable.long)
}


loadTopHitsQtl <- function(dirToSaveResultsQtlCombine,qtlDir, allQtlFiles,cytosToExclude=NULL, reExtractTopHits=F){
  ###Extract the top hits and create a single df for later processing
  
  if (reExtractTopHits){
    totalQtlRes = NULL 
    for (fileSel in allQtlFiles){
      print(fileSel)
      qtlRes = fread(file.path(qtlDir,fileSel), nrows = 100000)
      print(head(qtlRes))
      signifIndices = which(qtlRes$pvalue<1e-3)
      if (length(signifIndices)>0){
        qtlResSel = qtlRes[signifIndices,,drop=F]
        if (is.null(totalQtlRes)){
          totalQtlRes = qtlResSel
        }else{
          totalQtlRes = rbind.data.frame(totalQtlRes,qtlResSel)
          totalQtlRes = totalQtlRes[order(totalQtlRes$pvalue,decreasing = F),]
        }
      }
      print(head(totalQtlRes))
    }
    
    #Per cytokine define which SNPs to potentially annotate
    totalQtlRes$cytokine = totalQtlRes$gene
    snpsToAnnotate = defineSnpsToPotentiallyAnnotate(convertQtlToCorrectFormatForPlot(totalQtlRes),20,5e-8)
    fwrite(snpsToAnnotate, file = file.path(dirToSaveResultsQtlCombine,'snpsToAnnotate.txt'))
    
    print('SNPs to annotate')
    print(snpsToAnnotate)
    print(head(totalQtlRes$snps))
    totalQtlRes$annotate = 0
    totalQtlRes$annotate[which(totalQtlRes$snps %in% snpsToAnnotate$SNP)] = 1
    ###In the end I do not use this annotation, instead I use the one defined globally
    fwrite(x = totalQtlRes, file = file.path(dirToSaveResultsQtlCombine,paste0('combinedTopHits_',subFolderForSave,'.csv')))
  }else{
    totalQtlRes = fread(file.path(dirToSaveResultsQtlCombine,paste0('combinedTopHits_',subFolderForSave,'.csv')))
  }
  print(head(totalQtlRes))
  ##Exclude some cytokines from the plots
  if (!is.null(cytosToExclude)){
  rowsToKeep = grep(pattern = paste0('(',paste(cytosToExclude,collapse='|'),')'),totalQtlRes$cytokine, invert = T)
  print(head(rowsToKeep))
  totalQtlRes = totalQtlRes[rowsToKeep,]
  }
  return(totalQtlRes)
}

loadDosagesSnpList <- function(selSnps,folderGenotypes,chromsToChecks=c(1:22)){
  totalDosages = NULL
  for (chrSel in chromsToChecks){
    print(chrSel)
    rdsOutput=file.path(folderGenotypes,'RDS','mafHweOutlierFilteredGrCh38')
    saveFileGeneticsChrom = file.path(rdsOutput,paste0("GRCh38_300BCG_chr_",chrSel,"_annotated_MAF0.1_HWE1e-5_noOutliers.recode.tags_dosages.RDS"))
    genotypesDosages = readRDS(saveFileGeneticsChrom)
    intersectingSnps = intersect(selSnps,rownames(genotypesDosages))
    if (length(intersectingSnps)>0){
      genotypesDosages.sel = genotypesDosages[intersectingSnps,,drop=F]
    }else{
      next
    }
    if (is.null(totalDosages)){
      totalDosages = genotypesDosages.sel
    }else{
      totalDosages = rbind.data.frame(totalDosages,genotypesDosages.sel)
    }
  }
  return(totalDosages)
}


selectNonCorrelatedVariables <- function(dataframe.sel,sqCorCutoff=0.8){
  ##VARS SHOULD BE ORDERED FROM MOST TO LEAST SIGNIFICANT
  correlationsVars = cor(dataframe.sel, use='pairwise.complete.obs')^2
  varsToRemove = c()
  for (varIndex in c(1:nrow(correlationsVars))){
    varSel = rownames(correlationsVars)[[varIndex]]
    if (varSel %in% varsToRemove | varIndex==nrow(correlationsVars)){next}
    allCors = correlationsVars[c((varIndex+1):nrow(correlationsVars)),varSel]
    toRemoveExtra = names(allCors)[allCors>=sqCorCutoff]
    varsToRemove = c(varsToRemove,toRemoveExtra)
  }
  varsToKeep = setdiff(colnames(dataframe.sel),varsToRemove)
  return(varsToKeep)
}


compareTwoModels <- function(dataFrameSel, formula1, formula2){
  #Compare the fit of two linear models using ANOVA
  #Output: pvalue and explained variance
  
  #Remove any missing
  allIndependentVar = labels(terms(formula(formula2)))
  dataFrameSel.nomissing = dataFrameSel[which(rowSums(is.na(dataFrameSel[,allIndependentVar]))==0),]
  
  modelA = lm(formula = paste0(str_split(formula1, pattern='~')[[1]][[1]], ' ~ 1') , data = dataFrameSel.nomissing ) #basic formula with just the outcome variable
  modelB = lm(formula = formula1, data = dataFrameSel.nomissing) #everything except for the extra factors
  modelC = lm(formula = formula2, data = dataFrameSel.nomissing) #all factors
  ##Compare fit of the two models
  anovaRes = anova(modelA, modelB, modelC)
  
  totalVar = anovaRes[1,'RSS']
  B.perc = 100*anovaRes[2,'Sum of Sq']/totalVar
  C.perc = 100*anovaRes[3,'Sum of Sq']/totalVar
  C.pval = anovaRes[3,'Pr(>F)']
  return(list('pval'=C.pval,'varExpl'=C.perc))
}