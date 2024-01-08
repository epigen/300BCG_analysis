##This script prepares data to run the QTL enrichment using the PASCAL tool and subsmits the job to be run on a slurm based cluster
##The matching .batch file is "300BCG_QTLEnrichment.batch", which runs the actual PASCAL tool.
## 

#The original QTL analysis was performed on GrCh38, however pascal only runs on Hg19
#This is not an issue, since the rs-ids follow their own versioning independent of
#genome version (mostly), only the position of the SNP changes.
#We therefore ran pascal on GrCh37/HG19 
#We used custom gene sets. Specifically those used in EnrichR.

library(R.utils)
library(data.table)
source('processingFxnsGeneral.R')
source('plottingFxnsGeneral.R')
source('300bcg_generalFxns.R')

###Get the variables we need
library(config)
Sys.setenv(R_CONFIG_ACTIVE = "300bcg")
config <- config::get()
for (x in names(config)){eval(parse(text=paste0(x,"='",config[[x]],"'")))}

#Get the system arguments
args = base::commandArgs(trailingOnly=TRUE)

###Examples on how to run the script
# Rscript 300BCG_QTLEnrichment.R cytoData.good.log2 V1 on
# Rscript 300BCG_QTLEnrichment.R cytoData.good.fc V3 on
# Rscript 300BCG_QTLEnrichment.R cytoData.good.fc V2 on
# Rscript 300BCG_QTLEnrichment.R circMarkerData.log2 V1 on

###Example arguments
#args = c('cytoData.good.log2','V1','on')
#args = c('cytoData.good.fc','V3','on')

print(paste0('arguments=',args))
print(length(args))

#Load cytokine and metadata
allRelevantDataList = loadLukasMetaDataPbmc(allMetaDataPath)
list2env(allRelevantDataList, .GlobalEnv)

dataNameSel = args[[1]]
nameToAddToFolder = dataNameSel
selData = eval(parse(text = dataNameSel))
selVariables = colnames(selData) 

print('these variables are selected')
print(selVariables)

mainTpSel = args[[2]]
runPathwayEnrichment = args[[3]]

#In the end I am always running the analysis per chromosome to save memory
splitSetting = 'perChrom'

tpSel = mainTpSel
subFolderForSave = paste('quantitative',dataNameSel,paste(tpSel,collapse='.'),sep='__')

#Define where the QTL results are located
folderWithQtlFiles = file.path(genetics_save_dir,'QTLresults',subFolderForSave,splitSetting)
allQtlFiles = list.files(folderWithQtlFiles, pattern = ".*_withInfo.cqtls")

#Check if all files are there - if not - exit
selData = allRelevantDataList[[dataNameSel]]
allVariables = colnames(selData)
combisDf = expand.grid(allVariables,c(1:22))
allCombis = sort(paste0(combisDf[,1],"_chr",combisDf[,2],"_withInfo.cqtls"))
stopifnot(length(setdiff(allCombis,allQtlFiles))==0)

saveFolderPascalFiles = file.path(data_dir_300BCG,'genetics',"pascalInputFiles",subFolderForSave)
mkdirs(saveFolderPascalFiles)

geneSetFiles =  c("KEGG_2019_Human","GO_Biological_Process_2018")
pascalMethods = c('max')

#To run this on SLURM create a file with all combinations we want to run
#And later run these in a job array
#With combinations I mean: geneSet, max/sum, cytokines
if (runPathwayEnrichment=='on'){
  setMethodVariableCombis = expand.grid(list('geneSets'=geneSetFiles, 'pascalMethods'=pascalMethods,'allVariables'=allVariables))
  ncombis = nrow(setMethodVariableCombis)
  fileWithCombis = file.path(saveFolderPascalFiles,'allCombis.csv')
  fwrite(x = setMethodVariableCombis, file = fileWithCombis, sep=' ', col.names = F)
  addToName = paste0('300BCG__',subFolderForSave)
}else{
  setMethodVariableCombis = expand.grid(list('pascalMethods'=pascalMethods,'allVariables'=allVariables))
  setMethodVariableCombis = cbind(geneSets='none',setMethodVariableCombis)
  ncombis = nrow(setMethodVariableCombis)
  fileWithCombis = file.path(saveFolderPascalFiles,'allCombis_noEnrichment.csv')
  fwrite(x = setMethodVariableCombis, file = fileWithCombis, sep=' ', col.names = F)
  addToName = paste0('300BCG__',subFolderForSave,"_noEnrich")
}


#Now create the input tdf files that are needed to actually run PASCAL for all combinations
varIndex = 0
for (cytokineSel in allVariables){
  varIndex = varIndex+1
  print(paste(varIndex, " out of ", length(allVariables)))
  print(cytokineSel)
  #Load the QTL results per chromosome
  if (splitSetting=="perChrom"){
    extractedQtlData = NULL
    for (chrSel in c(1:22)){
      print(paste0('loading data for chromosome ', chrSel))
      qtlFileSel = paste0(cytokineSel,"_chr",chrSel,"_withInfo.cqtls")
      qtlData = fread(file.path(folderWithQtlFiles,qtlFileSel))
      extractedQtlDataPart = qtlData[,c('snps','pvalue')]
      if (is.null(extractedQtlData)){
        extractedQtlData = extractedQtlDataPart
      }else{
        extractedQtlData = rbind.data.frame(extractedQtlData,extractedQtlDataPart)
      }
    }
  }else{
    ###This part is basically depricated, we do not really use it anymore
    qtlFileSel = paste0(cytokineSel,"_withInfo.cqtls")
    qtlData = fread(file.path(folderWithQtlFiles,qtlFileSel))
    extractedQtlData = qtlData[,c('snps','pvalue')]
  }
  
  #now order by p-value
  extractedQtlData = extractedQtlData[order(extractedQtlData$pvalue, decreasing=F),]
  
  #Keep just those with an rs-id (and not those based on chrom-pos, since these do not translate to Hg19 easily)
  nonRs = which(!grepl(pattern='^rs',extractedQtlData$snps))
  withRs = which(grepl(pattern='^rs',extractedQtlData$snps))
  extractedQtlData = extractedQtlData[withRs,]
  #Write the files
  inputFilePascal = file.path(saveFolderPascalFiles,paste0(cytokineSel,'_pascal.tdf'))
  fwrite(extractedQtlData, file = inputFilePascal, col.names = F, sep = '\t')
}

#Get current filename and path (not very pretty code, maybe improve at some point)
#This way I can copy this code to other files with a matching batch file and it will still work
initial.options <- commandArgs(trailingOnly = FALSE)
file.arg.name <- "--file="
script.name <- sub(file.arg.name, "", initial.options[grep(file.arg.name, initial.options)])
script.basename <- system('pwd', intern = T)
print(script.basename)
print(script.name)
errorDir = file.path(data_dir_300BCG,'stdOutErr')
mkdirs(errorDir)

#WE USE HG19 here in PASCAL because that is what pascal was made for. However, it should not be a problem to use the
#GrCh38 QTL results, since rs-ids are stable between genome versions

if (runPathwayEnrichment=='on'){
  #Create a job array to run the pascal jobs,
  #here we convert per SNP p-values to per gene p-values AND
  #also get pathway enrichment
  command = paste0("sbatch", 
                   " --job-name='",subFolderForSave,"_pathways'",
                   " --error=",  "'",file.path(errorDir,subFolderForSave),"'",
                   ' --array=1-',ncombis,'%300 ',
                   " --export=fileWithCombis='",fileWithCombis,
                   "',runPathwayEnrichment='", runPathwayEnrichment,
                   "',outsuffix='", addToName,
                   "',inputFolderPascal='", saveFolderPascalFiles,"' ",
                   file.path(script.basename, gsub('\\.R$','.batch',script.name)))
}else{
  #here we just convert per SNP p-values to per gene p-values (if we want to GSEA ourselves)
  command = paste0("sbatch", 
                   " --job-name='",subFolderForSave,"_pathwaysNoEnrichment'",
                   " --error=",  "'",file.path(errorDir,subFolderForSave),"'",
                   ' --array=1-',ncombis,'%300 ',
                   " --export=fileWithCombis='",fileWithCombis,
                   "',outsuffix='", addToName,
                   "',runPathwayEnrichment='", runPathwayEnrichment,
                   "',inputFolderPascal='", saveFolderPascalFiles,"' ",
                   file.path(script.basename, gsub('\\.R$','.batch',script.name)))
}

#Actually run the command
print(command)
system(command)