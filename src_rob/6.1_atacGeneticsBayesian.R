library("R.utils")
library('openxlsx')
library(ggplot2)
library(stringr)
library(data.table)
library(dplyr)
library(stats)
library(BayesRRcpp)
library(R.utils)

###Get the variables we need
library(config)
Sys.setenv(R_CONFIG_ACTIVE = "300bcg")
config <- config::get()
for (x in names(config)){eval(parse(text=paste0(x,"='",config[[x]],"'")))}

source('300bcg_generalFxns.R')

regressOutCovariatesFromData <- function(inputData, covariateData, correctionFactors, visitSelY, visitSelCovariates='V1'){
  #Regress out the covariates from ATAC data
  #inputData - the atac data that is already normalized beforehand (these should contain just the regions I want to look at)
  #covariateData - data to be used or the covariate correction
  #correctionFactors - define the factors to correct the data for (these should be in "covariateData")
  #visitSel - visit to filter out
  
  allVars = colnames(inputData)
  covariateData$rn = rownames(covariateData)
  covariateData = covariateData[grep(paste0('.*\\_',visitSelCovariates,'$'),covariateData$rn),,drop=F]
  covariateData$rn = gsub(paste0(visitSelCovariates,'$'),visitSelY,covariateData$rn)
  
  inputData.resid = inputData
  inputData$rn = rownames(inputData)
  inputData = inputData[grep(paste0('.*\\_',visitSelY,'$'),rownames(inputData)),,drop=F]
  
  inputAndCovariates = base::merge(inputData,covariateData[,c(setdiff( colnames(covariateData),colnames(inputData)),'rn'), drop=F],by='rn')
  rownames(inputAndCovariates) = inputAndCovariates$rn
  
  inputData.resid[] = NA
  indexSel = 0
  nrOfRegions = length(allVars)
  for (varSel in allVars){
    indexSel = indexSel + 1
    print(paste0(indexSel, " out of ",nrOfRegions))
    formulaSel = paste0(varSel, ' ~ ', paste(correctionFactors, collapse=' + '))
    lm.res = lm(formula = as.formula(formulaSel), data = inputAndCovariates)
    residuals.sel = resid(lm.res)
    inputData.resid[names(residuals.sel),varSel] = residuals.sel
  }
  #Make sure they are ordered from most to least signif
  inputData.resid = inputData.resid[,allVars,drop=F]
  inputData.resid$rn = rownames(inputData.resid)
  #Now I already select for visit 1 (d0) sooner
  return(inputData.resid)
}

mapCytokineNamesForPlot <- function(oldLabs, color.by){
  if (color.by=='stimulus'){
    cytokinesSel = str_split_fixed(oldLabs, pattern = '_', n = 5)[,4]
    mappedCytokines = unlist(lapply(cytokinesSel, function(x) mappingListCyto[[x]]))
    stimsShortMapped = mapLongCytoToShortCytoName(str_split_fixed(oldLabs, pattern = '_', n = 5)[,1])
    newLabs = paste0( mappedCytokines, ' (', stimsShortMapped, ')')
  }else {
    newLabs = unlist(lapply(oldLabs, function(x) mappingListCircMed[[x]]))
  }
  return(newLabs)
}

args = base::commandArgs(trailingOnly=TRUE)

rdsOutput=file.path(genetics_save_dir,'RDS','mafHweOutlierFilteredGrCh38')
folderWithVcf = file.path(genetics_save_dir,"vcf","mafHweOutlierFilteredGrCh38")
folderForBedOutput = file.path(genetics_save_dir,"plink","mafHweOutlierFilteredGrCh38")

#Examples
#sbatch --job-name='bayesCytoV1_TTFT_20' --array=1-28%28 300bcg_atacGeneticsBayesian.batch cytoData.good.log2 V1 14 1000000 500000 T T F T 20
#args = c("cytoData.good.log2", "V1",'1','1000000','900000','T','T','F','T','50','0','0','28')

if (length(args)<2){
  
  for (seedSel in c("1","2","3","4","5")){
    nrOfVars="50"
    regressOut = "T"
    geneticsPCutoff = "0"
    atacPCutoff = "0"
    system(paste0("sbatch --job-name='bayesCytoV1_",regressOut,regressOut,"FT_",nrOfVars,"_",seedSel,"' --array=1-28%10 300bcg_atacGeneticsBayesian.batch cytoData.good.log2 V1 ",
                  seedSel," 1000000 900000 ",regressOut," ",regressOut," F T ",nrOfVars," ",atacPCutoff," ",geneticsPCutoff))
    
    system(paste0("sbatch --job-name='bayesCircV1_",regressOut,regressOut,"FT_",nrOfVars,"_",seedSel,"' --array=1-73%40 300bcg_atacGeneticsBayesian.batch circMarkerData.log2 V1 ",
                  seedSel," 1000000 900000 ",regressOut," ",regressOut," F T ",nrOfVars," ",atacPCutoff," ",geneticsPCutoff))
    
    system(paste0("sbatch --job-name='bayesCytoV3_",regressOut,regressOut,"FT_",nrOfVars,"_",seedSel,"' --array=1-28%10 300bcg_atacGeneticsBayesian.batch cytoData.good.fcCorr V3 ",
                  seedSel," 1000000 900000 ",regressOut," ",regressOut," F T ",nrOfVars," ",atacPCutoff," ",geneticsPCutoff))
    
    system(paste0("sbatch --job-name='bayesCytoV2_",regressOut,regressOut,"FT_",nrOfVars,"_",seedSel,"' --array=1-28%10 300bcg_atacGeneticsBayesian.batch cytoData.good.fcCorr V2 ",
                  seedSel," 1000000 900000 ",regressOut," ",regressOut," F T ",nrOfVars," ",atacPCutoff," ",geneticsPCutoff))
  }
}else if (length(args)==13){

  #arguments examples
  #args = c("cytoData.good.fcCorr","V3",'1','1000000','900000','T','T','F','T','50','0','0','27')
  #args = c("circMarkerData.log2","V1",'1','1000000','900000','T','T','F','T','50','0','0','27')
  #args = c("cytoData.good.log2", "V1",'1','1000000','900000','T','T','F','T','50','0','0','28')
  
  dataTypeSel = args[[1]]
  timepointSel = args[[2]]
  seedSel = as.numeric(args[[3]])
  maxIter = as.numeric(args[[4]])
  burnIn = as.numeric(args[[5]])
  regressOutOfAtac = as.logical(args[[6]])
  regressOutOfCells = as.logical(args[[7]])
  regressOutOfY = as.logical(args[[8]])
  filterTopHits = as.logical(args[[9]])
  topNrHits = as.numeric(args[[10]])
  atacPCutoff = 10^-as.numeric(args[[11]])
  geneticsPCutoff = 10^-as.numeric(args[[12]])
  
  factorIndex = as.numeric(args[[13]])
  
  print('args=')
  print(args)
  
  ##Folder to save the results
  subFolderSel = 'run1' #use the date to discrimate runs (not deterministic)
  folderToSaveBayesFiles = file.path(data_dir_300BCG,'genetics','atacGeneticsCompared','bayesRR',subFolderSel,paste(dataTypeSel, timepointSel,sep='_'))
  mkdirs(folderToSaveBayesFiles)
  
  folderToSaveManualFiles = file.path(data_dir_300BCG, 'genetics', 'atacGeneticsCompared','linearModelPercVariance',subFolderSel,paste(dataTypeSel, timepointSel,sep='_'))
  mkdirs(folderToSaveManualFiles)
  
  ##Load all metadata and data (the outcome factor)
  allRelevantDataList = loadLukasMetaDataPbmc(allMetaDataPath)
  list2env(allRelevantDataList, .GlobalEnv)
  
  #Extract the factors of interest
  #In case of the fold-changes, we want to regress out the effects of changes in cell counts and the 
  #changes in time o visit beforehand. This data was prepared already
  selVar = colnames(allRelevantDataList[[dataTypeSel]])[[factorIndex]]
  selY = allRelevantData[,selVar,drop=F]
  selY = selY[which(grepl(pattern = paste0('.*\\_',timepointSel,'$'), x = rownames(selY))),,drop=F]
  
  ##Define host factors to estmate for the host and environmental factors
  selectedHostEnvFactors = c("SEX", "AGE", "BMI", "oralContraceptivesIncludingMen", "VISIT_TIME_REAL", "alcoholInLast24h")
  seasonVars = c("VISIT_DATE_2PI_COS","VISIT_DATE_2PI_SIN")
  if (grepl(pattern =  ".*\\.(fcCorr|fc)$", x = dataTypeSel)){
    #In the case of fold-changes we add season
    selectedHostEnvFactors = c(selectedHostEnvFactors,seasonVars)
    
    #In that case we also remove evening individuals
    toRemoveEve = rownames(allRelevantData)[allRelevantData$IC_EVENING]
  }else{
    toRemoveEve = NULL
  }
  #Convert to numeric
  selectedHostEnv.df = allRelevantData[,selectedHostEnvFactors]
  selectedHostEnv.df$SEX = as.numeric(as.factor(selectedHostEnv.df$SEX))
  selectedHostEnv.df$alcoholInLast24h = as.numeric(selectedHostEnv.df$alcoholInLast24h)
  
  #Take the timepoint of interest
  selectedHostEnv.df = selectedHostEnv.df[which(grepl(pattern = paste0('.*\\_','V1','$'), x = rownames(selectedHostEnv.df))),]
  
  colnames(selectedHostEnv.df) = chartr(':/','..',colnames(selectedHostEnv.df))
  
  ##Define which cell counts to use for the estimation of variance of the cell counts
  correctionFactorList = loadCorrectionFactorList()
  cellVarsSel_V1 = correctionFactorList[[gsub('\\.(fcCorr|fc)$','.log2',dataTypeSel)]]
  cellVarsSel_V1 = setdiff(cellVarsSel_V1,selectedHostEnvFactors)
  cellVars.df = allRelevantData[,cellVarsSel_V1]
  
  ####Regress out covariates from Y
  if (regressOutOfY){
    ##Define correction
    correctionFactors = c(selectedHostEnvFactors, cellVarsSel_V1)
    selCovariateData = allRelevantData[,correctionFactors]#allRelevantData[,correctionFactors,drop=F]
    selCovariateData$SEX = as.numeric(as.factor(selCovariateData$SEX))
    selCovariateData$alcoholInLast24h = as.numeric(selCovariateData$alcoholInLast24h)
    selY = regressOutCovariatesFromData(selY, selCovariateData, correctionFactors, visitSelY=timepointSel, visitSelCovariates='V1')
    selY$rn <- NULL
  }
  
  
  ##Load files to map stimuli and cytokines to prettier names
  mapCytokineNamesFile = file.path('StimulusToPrettyNames.xlsx')
  sheetName = 'cytokine'
  loadedData = openxlsx::read.xlsx(xlsxFile = mapCytokineNamesFile, sheet = sheetName, colNames=F)
  mappingListCyto = createMappingList(loadedData)
  
  sheetName = 'circMed'
  loadedData = openxlsx::read.xlsx(xlsxFile = mapCytokineNamesFile, sheet = sheetName, colNames=F)
  mappingListCircMed = createMappingList(loadedData)
  
  #########                                               #########
  #######                                                   #######
  #####      Load and select genetics and epigenetics         #####
  ###                                                           ###
  #                                                               #
  
  #####Select top hits from genetics and epigenetics
  #I remove the regions and SNPs with R2>.8 from the top hits, this way we get more or less independent regions and SNPs
  #In the end it is mainly or only SNPs that are removed, since ATAC regions generally do not have such strong correlatins
  removeCorrelatedSnpsAndRegions=T
  if (filterTopHits){
    #########Load epigenetic-cytokine associations and select top regions
    ##### (1) ATAC-seq #####
    ##Load normalized epigenetics 
    if (startsWith(x = dataTypeSel, prefix = 'circ')){
      covariateAddToGrep = 'bmi\\.oralContra'
    }else{
      covariateAddToGrep = ''
    }
    
    fullAtac.long = loadCombinedAtacAssociationResults(allRelevantDataList, dataTypeSel,timepointSel, 
                                                       folderWithBaselineAtacResults=folderWithBaselineAtacResults, 
                                                       folderWithTrainingAtacResults=folderWithTrainingAtacResults,
                                                       mapFactorNames=F, mapPeakNames=F, covariateAddToGrep = covariateAddToGrep)
    
    ########Select the regions
    selVar.short = gsub('(\\.fc|\\.fcCorr)','',selVar)
    fullAtac.long.sel = fullAtac.long[which(fullAtac.long[,definitionsForSupplement[['factor']]] == selVar.short),]
    fullAtac.long.sel = fullAtac.long.sel[order(fullAtac.long.sel[,definitionsForSupplement[['pval']]]),]
    #Now also filter the top regions for p-value
    fullAtac.long.sel = fullAtac.long.sel[which(fullAtac.long.sel[definitionsForSupplement[['pval']]]<=atacPCutoff),]
    
    #here I select extra top regions, since we still want to remove potentially correlated ones after
    topRegionsSel.extended = fullAtac.long.sel[c(1:(topNrHits*20)),definitionsForSupplement[['region']]]
    topRegionsSel.extended = topRegionsSel.extended[which(!is.na(topRegionsSel.extended))]
    if (removeCorrelatedSnpsAndRegions){
      #Remove the correlated SNPs and regions from the top hits, but make sure that we still have the correct number of top hits
      ##Load ATACseq data --> we need this to calculate the correlations
      normalizedAtac = loadCountsPbmcLukas(fileWithAtacData)
      normalizedAtac = as.data.frame(t(normalizedAtac))
      ###Now, remove correlated ones
      nonCorrelatedRegions = selectNonCorrelatedVariables(dataframe.sel = normalizedAtac[,topRegionsSel.extended,drop=F], sqCorCutoff = .5)
      topRegionsSel = nonCorrelatedRegions[c(1:topNrHits)]
    }else{
      topRegionsSel = topRegionsSel.extended[c(1:topNrHits)]
    }
    topRegionsSel = topRegionsSel[!is.na(topRegionsSel)]
    
    ##### (2) Genetics #####
    ##########Select top SNPs
    subFolderForSave = paste('quantitative',gsub('.fcCorr','.fc',dataTypeSel),paste(timepointSel,collapse='.'),sep='__')
    qtlDir  = file.path(genetics_save_dir,'QTLresults',subFolderForSave,'perChrom')
    allQtlFiles = list.files(qtlDir, pattern = ".*_withInfo.cqtls")
    dirToSaveResultsQtlCombine = file.path(dirname(qtlDir),'combinedResults')
    totalQtlRes = loadTopHitsQtl(dirToSaveResultsQtlCombine,qtlDir, allQtlFiles,cytosToExclude=cytosToExclude, reExtractTopHits=F)
    totalQtlRes[['gene']] = gsub('.fc$','',totalQtlRes[['gene']])
    totalQtlRes.sel = totalQtlRes[totalQtlRes[['gene']]==selVar.short,]
    totalQtlRes.sel = totalQtlRes.sel[order(totalQtlRes.sel[,'pvalue']),]
    
    ##Filter these for a certain pvalue
    totalQtlRes.sel = totalQtlRes.sel[which(totalQtlRes.sel$pvalue<=geneticsPCutoff),,drop=F]
    
    if (removeCorrelatedSnpsAndRegions){
      #Load the genetics data
      #Keep only the significant snps above some threshold
      selSnps = unlist(totalQtlRes.sel[,'snps'])
      totalDosages = NULL
      for (chrSel in c(1:22)){
        print(chrSel)
        saveFileGeneticsChrom = file.path(rdsOutput,paste0("GRCh38_300BCG_chr_",chrSel,"_annotated_MAF0.1_HWE1e-5_noOutliers.recode.tags_dosages.RDS"))
        genotypesDosages = readRDS(saveFileGeneticsChrom)
        genotypesDosages = genotypesDosages[intersect(selSnps,rownames(genotypesDosages)),]
        if (is.null(totalDosages)){
          totalDosages = genotypesDosages
        }else{
          totalDosages = rbind.data.frame(totalDosages,genotypesDosages)
        }
      }
      totalDosages = t(totalDosages)
      rownames(totalDosages) = paste0(rownames(totalDosages),'_V1')
      totalDosages = totalDosages[,selSnps]
      
      nonCorrelatedRegions = selectNonCorrelatedVariables(dataframe.sel = totalDosages, sqCorCutoff = .5)
      
      topSnpsSel = nonCorrelatedRegions[c(1:topNrHits)]
      topSnpsSel = topSnpsSel[!is.na(topSnpsSel)]
      
    }else{
      topSnpsSel = unlist(totalQtlRes.sel[c(1:topNrHits),'snps'])
    }
  }else{
    #If we do not select the top X regions and SNPs, we need to get rid of correlated SNPs, because otherwise there are simply
    #too many for the algorithm to handle.
    
    #clump genetics data to only keep independent snps
    #https://privefl.github.io/bigsnpr/articles/pruning-vs-clumping.html
    library(bigsnpr)
    # Convert vcf to bed

    mkdirs(folderForBedOutput)
    oldWd = getwd()
    setwd(folderForBedOutput)
    prunedSnpList = c()
    for (chrSel in c(1:22)){
      print(chrSel)
      basicFileName = paste0('GRCh38_300BCG_chr_',chrSel,'_annotated_MAF0.1_HWE1e-5_noOutliers.recode.tags')
      rdsPath = file.path(folderForBedOutput, paste0(basicFileName,'.rds'))
      if (!file.exists(file.path(folderForBedOutput, paste0(basicFileName,'.rds')))){
        vcfFileSel = file.path(folderWithVcf,paste0(basicFileName, '.vcf.gz'))
        commandSel= paste0('~/local/plink --vcf ',vcfFileSel,' --make-bed --out ', paste0(basicFileName))
        system(commandSel)
        #Read the bed file
        rdsPath = snp_readBed(file.path(folderForBedOutput, paste0(basicFileName,'.bed')))
      }
      
      snpData = readRDS(rdsPath)
      indicesToKeep = snp_clumping(
        G=snpData$genotypes,
        infos.chr = snpData$map$chromosome,
        thr.r2 = 0.2,
        ncores = 1)
      
      snpData.map.pruned = snpData$map[indicesToKeep,]
      nonMissingRsIdIndices = which(snpData.map.pruned[,'marker.ID']!='.')
      missingRsIdIndices = which(snpData.map.pruned[,'marker.ID']=='.')
      rsSnpsToKeep = snpData.map.pruned[nonMissingRsIdIndices,'marker.ID']
      #For the non-rs ids we need to include both ways (first allele1 or allele2) of defining the alternative and reference allele, since this
      #order does not match with the way it was loaded with from the vcf. If we include both, one of the two will always overlap
      nonRsSnpsToKeep = c(paste0("chr",snpData.map.pruned[missingRsIdIndices,'chromosome'],
                                 ":",
                                 snpData.map.pruned[missingRsIdIndices,'physical.pos'],
                                 "_",
                                 snpData.map.pruned[missingRsIdIndices,'allele1'],"/",
                                 snpData.map.pruned[missingRsIdIndices,'allele2']),  
                          paste0("chr",snpData.map.pruned[missingRsIdIndices,'chromosome'],
                                 ":",
                                 snpData.map.pruned[missingRsIdIndices,'physical.pos'],
                                 "_",
                                 snpData.map.pruned[missingRsIdIndices,'allele2'],"/",
                                 snpData.map.pruned[missingRsIdIndices,'allele1']))
      snpsToKeep = c(rsSnpsToKeep, nonRsSnpsToKeep)
      prunedSnpList = c(prunedSnpList, snpsToKeep)
    }
    setwd(oldWd)
  }
  
  ###Define which set of SNPs to use
  if (filterTopHits){
    selSnpsForAnalysis = topSnpsSel
  }else{
    selSnpsForAnalysis = prunedSnpList
  }
  
  #Load the genetics data
  totalDosages = NULL
  for (chrSel in c(1:22)){
    print(chrSel)
    rdsOutput=file.path(genetics_save_dir,'RDS','mafHweOutlierFilteredGrCh38')
    saveFileGeneticsChrom = file.path(rdsOutput,paste0("GRCh38_300BCG_chr_",chrSel,"_annotated_MAF0.1_HWE1e-5_noOutliers.recode.tags_dosages.RDS"))
    genotypesDosages = readRDS(saveFileGeneticsChrom)
    genotypesDosages = genotypesDosages[intersect(selSnpsForAnalysis,rownames(genotypesDosages)),,drop=F]
    if (is.null(totalDosages)){
      totalDosages = genotypesDosages
    }else{
      totalDosages = rbind.data.frame(totalDosages,genotypesDosages)
    }
  }
  totalDosages = t(totalDosages)
  rownames(totalDosages) = paste0(rownames(totalDosages),'_V1')
  colnames(totalDosages) = chartr(':/','..',colnames(totalDosages))
  
  #########                                               #########
  #######                                                   #######
  #####     Regress out covariates epigenetics&cell count     #####
  ###                                                           ###
  #                                                               #
  folderToSaveRegressedOutData = file.path(data_dir_300BCG, "atacSeq")
  mkdirs(folderToSaveRegressedOutData)
  fileToSaveRegressedOutData = file.path(folderToSaveRegressedOutData,'regressed_out_batch_corrected_normalized_log2_CPM_PBMC.RDS')
  
  if (regressOutOfAtac){
    #Redo the regression or load the saved one
    reRegressOutCovariatesAtac = F
    if (reRegressOutCovariatesAtac){
      ##Load ATACseq data
      normalizedAtac = loadCountsPbmcLukas(fileWithAtacData)
      normalizedAtac = as.data.frame(t(normalizedAtac))
      ##Regress out the covariates ATAC
      ##The covariates for correcting the ATAC should always be based on V1
      #So make sure to get those covariates, and not the fold-change/d90 covariates
      correctionFactorsV1 =  c(cellVarsSel_V1,selectedHostEnvFactors)#correctionFactorList[['cytoData.good.log2']] #This set of correction factors is what we want to use: we want to correct ATAC at baseline, plus sex, age and visit time
      atacExtraCorrections = c('TSS_ENRICHMENT')
      correctionFactors.atac = unique(c(cellVarsSel_V1,selectedHostEnvFactors,atacExtraCorrections))
      selCovariateData.atac = allRelevantData[,correctionFactors.atac]
      
      normalizedAtac.resid = regressOutCovariatesFromData(normalizedAtac, covariateData = selCovariateData.atac, 
                                                          correctionFactors = correctionFactors.atac, visitSelY='V1', visitSelCovariates='V1')
      saveRDS(object = normalizedAtac.resid, file = fileToSaveRegressedOutData)
    }else{
      normalizedAtac.resid = readRDS(file = fileToSaveRegressedOutData)
    }
    normalizedAtac.resid$rn <- NULL
  }else{
    normalizedAtac = loadCountsPbmcLukas(fileWithAtacData)
    normalizedAtac = as.data.frame(t(normalizedAtac))
    normalizedAtac.resid = normalizedAtac
  }
  
  ###Select top atac hits
  if (filterTopHits){
    normalizedAtac.resid = normalizedAtac.resid[,topRegionsSel,drop=F]
  }
  colnames(normalizedAtac.resid) = chartr(':/','..',colnames(normalizedAtac.resid))
  
  #select columns with cell counts
  if (regressOutOfCells){
    #Regress age, sex, and the other covariates out of the cell type data
    selCovariateData.cellTypes = allRelevantData[,selectedHostEnvFactors,drop=F]
    cellVars.df.resid = regressOutCovariatesFromData(cellVars.df, covariateData = selCovariateData.cellTypes, 
                                                     correctionFactors = selectedHostEnvFactors, 
                                                     visitSelY="V1", visitSelCovariates='V1')
    cellVars.df[,cellVarsSel_V1] = cellVars.df.resid[,cellVarsSel_V1] 
  }
  colnames(cellVars.df) = chartr(':/','..',colnames(cellVars.df))
  
  #########                                               #########
  #######                                                   #######
  #####            Merge and scale all the data               #####
  ###                                                           ###
  #                                                               #
  
  #Merge these and 
  #(1) make sure there are no NAs
  #(2) make sure they are all scaled
  #Common ids
  
  rownames(selY) = gsub('\\_V(2|3)$','_V1',x = rownames(selY))
  commonIds = Reduce(intersect, list(rownames(selY),
                                     rownames(selectedHostEnv.df),
                                     rownames(cellVars.df),
                                     rownames(totalDosages),
                                     rownames(normalizedAtac.resid)))
  
  #Remove evening for the fold-changes
  if (!is.null(toRemoveEve)){
    commonIds = setdiff(commonIds,toRemoveEve)
  }
  
  #Merge data
  if (regressOutOfY){
    mergedAllData = cbind.data.frame(selY[commonIds,,drop=F], 
                                     normalizedAtac.resid[commonIds,,drop=F],
                                     totalDosages[commonIds,,drop=F])
  }else{
    mergedAllData = cbind.data.frame(selY[commonIds,,drop=F],
                                     selectedHostEnv.df[commonIds,,drop=F], 
                                     cellVars.df[commonIds,,drop=F], 
                                     normalizedAtac.resid[commonIds,,drop=F],
                                     totalDosages[commonIds,,drop=F] 
    )
  }
  #Remove NA's
  toKeepNoNa = which(rowSums(is.na(mergedAllData))==0)
  mergedAllData = mergedAllData[toKeepNoNa,,drop=F]
  
  #Split into X and Y
  mergedYData= mergedAllData[,which(colnames(mergedAllData)==selVar),drop=F]
  mergedXData= mergedAllData[,which(colnames(mergedAllData)!=selVar),drop=F]
  
  #Scale X and Y
  Y <- scale(mergedYData)
  X <- scale(mergedXData)
  
  #Make sure the colnames do not contain special characters
  colnames(X) = chartr(':/','..',colnames(X))
  
  #########                                               #########
  #######                                                   #######
  #####            Also do the manual calculation             #####
  ###                                                           ###
  #                                                               #
  
  ##Also do the manual calculation
  #(1) First look at host factors
  selFactorsHost = setdiff(colnames(selectedHostEnv.df),seasonVars)
  hostFactorFormula = paste0(selVar, ' ~ ' , paste(selFactorsHost, collapse= ' + '))
  lm.res = lm(mergedAllData, formula = hostFactorFormula)
  anova.res = anova(lm.res)
  percHost = 100*sum(anova.res[selFactorsHost,'Sum Sq']/sum(anova.res[,'Sum Sq']))
  print(percHost)
  
  #(1b)If we are looking at fold-changes, we also want season
  if (grepl(pattern =  ".*\\.(fcCorr|fc)$", x = dataTypeSel)){
    hostFactorFormula2 = paste0(hostFactorFormula, ' + ', paste(seasonVars, collapse = ' + '))
    results1 = compareTwoModels(mergedAllData, hostFactorFormula, hostFactorFormula2)
    percSeason = results1[['varExpl']]
    print(percSeason)
  }else{
    hostFactorFormula2 = hostFactorFormula
  }
  
  #(2) Cell counts
  formulaCellCounts = paste0(hostFactorFormula2, ' + ', paste(colnames(cellVars.df), collapse = ' + '))
  results2 = compareTwoModels(mergedAllData, hostFactorFormula2, formulaCellCounts)
  percCells = results2[['varExpl']]
  print(percCells)
  
  #(3) Atac-seq over genetics
  formulaAddedGenetics = paste0(formulaCellCounts, ' + ', paste(colnames(totalDosages), collapse = ' + '))
  formulaAddedGeneticsAtac = paste0(formulaAddedGenetics, ' + ', paste(colnames(normalizedAtac.resid), collapse = ' + '))
  results3 = compareTwoModels(mergedAllData, formulaAddedGenetics, formulaAddedGeneticsAtac)
  percAtac = results3[['varExpl']]
  print(percAtac)
  
  #(4) Genetics over atac
  formulaAddedAtac = paste0(formulaCellCounts, ' + ', paste(colnames(normalizedAtac.resid), collapse = ' + '))
  formulaAddedGeneticsAtac = paste0(formulaAddedAtac, ' + ', paste(colnames(totalDosages), collapse = ' + '))
  results4 = compareTwoModels(mergedAllData, formulaAddedAtac, formulaAddedGeneticsAtac)
  percGenetics = results4[['varExpl']]
  print(percGenetics)
  
  #(5) Residuals
  lm.res = lm(mergedAllData, formula = formulaAddedGeneticsAtac)
  anova.res = anova(lm.res)
  percResiduals = 100*sum(anova.res['Residuals','Sum Sq']/sum(anova.res[,'Sum Sq']))
  print(percResiduals)
  
  ##What if I do it for the total
  percHost2 = 100*sum(anova.res[colnames(selectedHostEnv.df),'Sum Sq']/sum(anova.res[,'Sum Sq']))
  percCells2 = 100*sum(anova.res[colnames(cellVars.df),'Sum Sq']/sum(anova.res[,'Sum Sq']))
  percAtac2 = 100*sum(anova.res[colnames(normalizedAtac.resid),'Sum Sq']/sum(anova.res[,'Sum Sq']))
  percGenetics2 = 100*sum(anova.res[colnames(totalDosages),'Sum Sq']/sum(anova.res[,'Sum Sq']))
  percResiduals2 = 100*sum(anova.res['Residuals','Sum Sq']/sum(anova.res[,'Sum Sq']))
  if (grepl(pattern =  ".*\\.(fcCorr|fc)$", x = dataTypeSel)){
    percSeason2 = 100*sum(anova.res[seasonVars,'Sum Sq']/sum(anova.res[,'Sum Sq']))
  }
  
  print(percHost2)
  print(percCells2)
  print(percAtac2)
  print(percGenetics2)
  print(percResiduals2)
  
  if (grepl(pattern =  ".*\\.(fcCorr|fc)$", x = dataTypeSel)){
    percExplainedAnova = rbind(c(percHost,percSeason,percCells,percAtac,percGenetics,percResiduals),
                               c(percHost2,percSeason2,percCells2,percAtac2,percGenetics2,percResiduals2))
    colnames(percExplainedAnova) = c('host','season','cellCounts','atac','genetics','residuals')
  }else{
    percExplainedAnova = rbind(c(percHost,percCells,percAtac,percGenetics,percResiduals),
                               c(percHost2,percCells2,percAtac2,percGenetics2,percResiduals2))
    colnames(percExplainedAnova) = c('host','cellCounts','atac','genetics','residuals')
  }
  
  rownames(percExplainedAnova) = c('manualOrder','sequential')
  fwrite(percExplainedAnova, file = file.path(folderToSaveManualFiles,paste0(selVar,"_",paste(args, collapse='_'),".csv")), row.names = T)
  
  
  #########                                               #########
  #######                                                   #######
  #####            Do the actual basyesian calc               #####
  ###                                                           ###
  #                                                               #
  
  for (geneticsInclusion in c('withGenetics','withoutGenetics')){
    if (regressOutOfY){
      G=2
      GroupAssignment = c(rep(0,ncol(normalizedAtac.resid)),
                          rep(1,ncol(totalDosages)))
    }else{
      #Separate groups for the host factor cocariates
      selectedHostEnvFactorsNoSeason = setdiff(selectedHostEnvFactors,seasonVars)
      GroupAssignmentPt1 = c(rep(0,length(selectedHostEnvFactorsNoSeason)))
      
      #In the case of fold-changes add season
      names(GroupAssignmentPt1) = selectedHostEnvFactorsNoSeason
      if (grepl(pattern =  ".*\\.(fcCorr|fc)$", x = dataTypeSel)){
        GroupAssignmentPt1b = c(rep(1,length(seasonVars)))
        names(GroupAssignmentPt1b) = seasonVars
        GroupAssignmentPt1 = c(GroupAssignmentPt1,GroupAssignmentPt1b)
      }
      
      GroupAssignmentPt2 = c(rep(max(GroupAssignmentPt1)+1,length(cellVarsSel_V1)))
      names(GroupAssignmentPt2) = cellVarsSel_V1

      GroupAssignment = c(GroupAssignmentPt1,GroupAssignmentPt2,
                          rep(max(GroupAssignmentPt2)+1,ncol(normalizedAtac.resid)),
                          rep(max(GroupAssignmentPt2)+2,ncol(totalDosages)))
      G=length(unique(GroupAssignment))
    }
    
    ##Remove the genetics part if that is indicated
    if (geneticsInclusion=='withoutGenetics'){
      nrToRemove = length(selSnpsForAnalysis)
      GroupAssignment.sel = GroupAssignment[c(1:(length(GroupAssignment)-nrToRemove))]
      G.sel =  length(unique(GroupAssignment.sel))
      X.sel =  X[,c(1:(length(GroupAssignment)-nrToRemove))]
    }else{
      X.sel=X
      G.sel = G
      GroupAssignment.sel = GroupAssignment
    }
    
    #P=0.5 #prior probability of a marker being excluded from the model
    sigma0=0.01# prior  variance of a zero mean gaussian prior over the mean mu NOT IMPLEMENTED
    v0E= 0.0001 # degrees of freedom over the inv scaled chi square prior over residuals variance
    s02E = 0.001 #scale of the inv scaled chi square prior over residuals variance
    v0G = 0.0001 #degrees of freedom of the inv bla bla prior over snp effects
    s02G = 0.001 # scale for the samecva 
    cva=matrix(rep(c(0.01,0.1, 1),G),nrow=G,byrow = T)  #components variance prior scales
    set.seed(seedSel)
    
    ## a parameter seed is present but has not effect as we rely on the R environment seed instead!!
    bayesROutputFile = file.path(folderToSaveBayesFiles,paste0(selVar,"_",paste(args, collapse='_'),'_',geneticsInclusion,".csv"))
    BayesRSamplerV2Groups(outputFile = bayesROutputFile,
                          seed = seedSel, 
                          max_iterations = maxIter, 
                          burn_in = burnIn,
                          thinning = 10,
                          X = X.sel, 
                          Y = Y,
                          sigma0 = sigma0,
                          v0E = v0E,s02E = s02E,
                          v0G = v0G,s02G = s02G,
                          cva = cva,
                          groups = G.sel,
                          gAssign = GroupAssignment.sel)
  }
}
