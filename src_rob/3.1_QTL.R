library(MatrixEQTL)
library(snpStats)
library(R.utils)
library(data.table)
library(stringr)
library(dplyr)

start_time <- Sys.time()

sliceQuantAndMetaDataForQtl <- function(dataXAxis, newMetaData, correctionFactors){
  #Slice the data
  cvrt.mat = as.matrix(data.frame(newMetaData[,correctionFactors,drop=F]))
  cvrtSliced = SlicedData$new( t(cvrt.mat) );
  cvrtSliced$ResliceCombined(2000);
  
  #Slice the data
  cyto.mat = as.matrix(dataXAxis)
  cytoSliced = SlicedData$new( t(cyto.mat) );
  cytoSliced$ResliceCombined(2000);
  return(list('cvrtSliced'=cvrtSliced,'cytoSliced'=cytoSliced))
}

###Get the variables we need
library(config)
Sys.setenv(R_CONFIG_ACTIVE = "300bcg")
config <- config::get()
for (x in names(config)){eval(parse(text=paste0(x,"='",config[[x]],"'")))}

source('processingFxnsGeneral.R')
source('plottingFxnsGeneral.R')
source('300bcg_generalFxns.R')

####EXAMPLES#####
#Rscript 300bcg_QTL.R cytoData.good.log2 V1
#Rscript 300bcg_QTL.R cytoData.good.fc V3
#Rscript 300bcg_QTL.R cytoData.good.fc V2
#Rscript 300bcg_QTL.R circMarkerData.log2 V1
#args = c('cytoData.good.log2',V1','1')
#args = c("circMarkerData.fc", "V3",'1')

#Get the system arguments
args = base::commandArgs(trailingOnly=TRUE)
#args = c('cytoData.good.fc','V3','1')
print(paste0('arguments=',args))
print(length(args))

#Load cytokine and metadata
#allMetaDataPath = "" #set path to data (complete_metadata.corrected.csv) as combined by Lukas 
allRelevantDataList = loadLukasMetaDataPbmc(allMetaDataPath)
list2env(allRelevantDataList, .GlobalEnv)

#Define the correctionfactors for the cell subtype part
correctionFactorList = loadCorrectionFactorList()

print("finished loading data")

################SETTINGS################
dataNameSel = args[[1]]
nameToAddToFolder = dataNameSel
selData = eval(parse(text = dataNameSel))
selVariables = colnames(selData) 
print('these variables are selected')
print(selVariables)
mainTpSel = args[[2]]

#Define the covariates
covariatesSel = correctionFactorList[[dataNameSel]]

#Define individuals to remove
if (endsWith(dataNameSel, '.fc')){
  #In the case of fold-changes at day 14 (V2) or day 90 (V3), we remove the indidivuals vaccinated in the evening
  toRemoveEve = rownames(allRelevantData)[allRelevantData$IC_EVENING]
}else{
  toRemoveEve = NULL
}

#Extract the data we need
################################
#NOTE: "allRelevantData" contains all the metadata/feature data we need, and we extract the variables from this
#Now create a data.frame with all combinations of variable and chromosome
variableDf = as.data.frame(selVariables)
variableDf$selVariables = as.character(variableDf$selVariables)
ncombis=nrow(variableDf)

#Either submit jobs or run the script
#If the number of input arguments is exactly 4, the QTL calculations will be run
#If the number is 3, jobs will be submitted with 4 arguments, with the added argument being an index. This index is used to define the quantitative variable of interest.
print('starting if-statement')
print(paste0('number of arguments is equal to ',length(args)))
if (length(args)<2 | length(args)>3){
  stop("wrong number input arguments")
}else if (length(args)==2){
  ##In case there are three arguments we submit the job array
  print('starting first part')
  #Run per chromosome to save memory
  perChromosome = T
  
  #Get current filename and path (not very pretty code, maybe improve at some point)
  initial.options <- commandArgs(trailingOnly = FALSE)
  file.arg.name <- "--file="
  script.name <- sub(file.arg.name, "", initial.options[grep(file.arg.name, initial.options)])
  script.basename <- system('pwd', intern = T)
  print(script.basename)
  print(script.name)
  errorDir = file.path(data_dir_300BCG,'stdOutErr')
  mkdirs(errorDir)
  
  origJobNamePart =  paste0(" --job-name='",nameToAddToFolder,"_QTLs'")
  origErrorFolderPart =  paste0(" --error=",  "'",file.path(errorDir,paste0(nameToAddToFolder)),".err'")
  
  basicCommand = paste0("sbatch", 
                        origJobNamePart,
                        origErrorFolderPart,
                        paste0(' --array=1-',ncombis,'%100'),
                        " --export=scriptName='",file.path(script.name),"'")
  fullScriptPath = file.path(gsub('\\.R$','.batch',script.name))
  command = paste0(basicCommand, ",dataNameSel='", dataNameSel,"'mainTpSel='", mainTpSel,"' ", fullScriptPath)
  print(command)
  system(command)
  
}else if(length(args)==3){
  #In case there are four arguments, this means that the job array is already constructed and we can run the actual job
  print('starting script')
  splitSetting = 'perChrom'
  print(paste0('splitSetting = ',splitSetting))
  
  #Check which specific one we are running
  variableNr = args[[3]]
  
  varSel = variableDf[[variableNr,1]]
  print(paste0('selected cytokine = ',varSel))
  for (chrSel in c(1:22)){
    #Go over all chromosomes
    #Note: I also tried to include chromosome as part of the array,
    #but since it takes some time for a job to start running, this meant
    #that getting all jobs done was taking much more time in that case
    print(paste0('selected chromesome = ',chrSel))
    
    #====SETTINGS====
    ##If we want to look at the absolute values or fold-changes we use the linear model
    tpSel = mainTpSel
    useModel = modelLINEAR
    
    #Keep just selected timepoint
    print(paste0("extract timepoint ",tpSel))
    print(head(allRelevantData))
    combinedData = allRelevantData[which(base::endsWith(rownames(allRelevantData),paste0('_',tpSel))),]
    subFolderForSave = paste('quantitative',dataNameSel,paste(tpSel,collapse='.'),sep='__')
    
    #Convert the covariate SEX into a numberic variable
    if ('SEX' %in% covariatesSel){
      combinedData$SEX[which(combinedData$SEX=='F')]=0
      combinedData$SEX[which(combinedData$SEX=='M')]=1
      combinedData$SEX = as.numeric(combinedData$SEX)
    }
    
    #Get the vars
    selCovariateData = combinedData[,covariatesSel]
    print('extracting numeric covariates')
    nums <- unlist(lapply(selCovariateData, is.numeric))
    selCovariateDataNum = as.matrix(selCovariateData[,nums])
    #Check that all variables have been converted to numeric
    stopifnot(length(setdiff(covariatesSel,colnames(selCovariateDataNum)))==0)
    selCovariateDataNum = selCovariateDataNum[rowSums(is.na(selCovariateDataNum))==0,]
    
    #Check that there are no linearly dependent variables
    print("define linearly dependent vars")
    linearlyIndependentToKeep = colnames(selCovariateDataNum)[qr(selCovariateDataNum)$pivot[seq_len(qr(selCovariateDataNum)$rank)]]
    linearlyDependentToRemove = setdiff(colnames(selCovariateDataNum),linearlyIndependentToKeep)
    print(paste0('removed ',paste(linearlyDependentToRemove,collapse = '; '), ' because of linear dependence'))
    linearlyIndependentToKeep = c(linearlyIndependentToKeep,setdiff(covariatesSel,colnames(selCovariateDataNum)))
    #Get them back to the right order 
    linearlyIndependentToKeep = intersect(covariatesSel,linearlyIndependentToKeep)
    #There should be not linearly dependent variables
    stopifnot(length(linearlyDependentToRemove)==0)
    
    #======Process genetics data======
    #load genetics data
    subFolderWithData='mafHweOutlierFilteredGrCh38'
    rdsOutput=file.path(genetics_save_dir,'RDS',subFolderWithData)
    
    saveFileGeneticsChrom = file.path(rdsOutput,paste0("GRCh38_300BCG_chr_",chrSel,"_annotated_MAF0.1_HWE1e-5_noOutliers.recode.tags_dosages.RDS"))
    genotypesDosages = readRDS(saveFileGeneticsChrom)
    variantInfo = readRDS(file.path(rdsOutput,paste0("chr_",chrSel,"_300bcg_variantInfoAddedChrPos.RDS")))
    
    #Define which variants to use
    #For now add '_V1' to the names of the individuals to match the structure of the metadata
    genotypesDosages_selVariants = genotypesDosages
    colnames(genotypesDosages_selVariants) = paste(colnames(genotypesDosages_selVariants),mainTpSel,sep='_')
    variantInfo_selVariants = variantInfo
    rm(genotypesDosages)
    
    ##The selected column is the variable we are interested in
    colSel = varSel
    print(paste0('extracting column ', colSel))
    #Slice the quantitative data and the metadata
    #We want to do this one cytokine at a time due to missing data sometimes being different for different variables
    #linearlyIndependentToKeep = the covariates
    quantData = combinedData[,colSel,drop=F]
    covariatesData = combinedData[,linearlyIndependentToKeep]
    print(head(quantData))
    print(head(covariatesData))
    
    #Now check for which samples we have complete data (no missing in any), and slice the matrices
    idsQuant = rownames(quantData)[which(base::rowSums(is.na(quantData))==0)]
    idsCov = rownames(covariatesData)[which(base::rowSums(is.na(covariatesData))==0)]
    idsSnps = colnames(genotypesDosages_selVariants)[which(base::colSums(is.na(genotypesDosages_selVariants))==0)]
    commonIds = Reduce(intersect, list(idsQuant,idsCov,idsSnps))
    
    #Remove evening sampes if they are there
    if (!is.null(toRemoveEve)){
      commonIds = setdiff(commonIds,toRemoveEve)
    }
    print(commonIds)
    
    #!!!CHECK IF OUTLIERS ARE REMOVED
    outliers = rownames(metaData)[which(rowSums(metaData[,c("EXCLUSION","SNP_OUTLIER")]!='')>0)]
    stopifnot(length(intersect(commonIds,outliers))==0)
    
    #Keep just the common ids between all data types
    print("keeping the common ids")
    quantData = quantData[commonIds,,drop=F]
    covariatesData = covariatesData[commonIds,]
    geneticsDataSel = genotypesDosages_selVariants[,commonIds]
    rm(genotypesDosages_selVariants)
    
    #Slice the SNPs as MagtrixEQTL wants
    print("slicing the SNPs")
    snpsSliced = SlicedData$new( as.matrix(geneticsDataSel) )
    snpsSliced$ResliceCombined(2000)
    
    sliceRes = sliceQuantAndMetaDataForQtl(quantData, covariatesData, linearlyIndependentToKeep)
    cvrtSliced = sliceRes[['cvrtSliced']]
    cytoSliced= sliceRes[['cytoSliced']]
    
    print(cvrtSliced)
    print(cytoSliced)
    print(snpsSliced)
    
    ##Output all data (i.e. pvOutputThreshold=1)
    output_file_name = ""
    pvOutputThreshold = 1
    
    ## Run the analysis
    matrix_cQTLS = Matrix_eQTL_engine(
      snps = snpsSliced,
      gene = cytoSliced,
      cvrt = cvrtSliced,
      output_file_name = output_file_name,
      pvOutputThreshold = pvOutputThreshold,
      useModel = useModel, 
      errorCovariance = numeric(), 
      verbose = FALSE,
      pvalue.hist = FALSE,
      min.pv.by.genesnp = FALSE,
      noFDRsaveMemory = FALSE);
    
    #Extract the relevant results
    fullResults <- matrix_cQTLS$all$eqtls
    
    #Save the results
    dirToSaveQtlOutput = file.path(genetics_save_dir,"QTLresults")
    folderToSaveResults = file.path(dirToSaveQtlOutput,subFolderForSave,splitSetting)
    addToFileNames = paste0('_chr',chrSel)
    mkdirs(folderToSaveResults)
    fullResults = fullResults[order(fullResults$pvalue, decreasing = F),]
    #merge this with the bim data
    #Add other info to the data. This only goes wrong for the SNPs without an RS-id, but I fix this in later scripts
    fullResWithInfo = dplyr::full_join(fullResults, variantInfo_selVariants,  by = c("snps" = "ID"))
    print(head(fullResWithInfo))
    fwrite(x = fullResWithInfo, file = file.path(folderToSaveResults,paste0(varSel,addToFileNames,'_withInfo.cqtls')))
  }
}else{
  stop('Tried to submit all jobs, something went wrong')
}

end_time <- Sys.time()
print(end_time - start_time)