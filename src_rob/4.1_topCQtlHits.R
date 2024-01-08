library(janitor)
library(R.utils)
library(data.table)
library(stringr)
library(dplyr)
library(ggplot2)
library(openxlsx)

createManhattanPlot <- function(gwasResults, snpsToAnnotate, pValCutoffMark1, pValCutoffMark2, colorByStims=T, 
                                textSize=12, dividingFactorNonSignifSnps = 15, colorOutline ='#937860', highlightSignificant=F,rasterizePlot=T){
  library(ggrepel)
  print('start normal manhattan')
  gwasResults$snp_trait = paste(gwasResults$SNP, gwasResults$trait, sep='__')
  snpsToAnnotate$snp_trait = paste(snpsToAnnotate$SNP, snpsToAnnotate$trait, sep='__')
  snp_trait_toAnnotate = snpsToAnnotate$snp_trait
  #=========First normal Manhattan=============
  # Prepare the dataset
  don <- gwasResults %>% 
    
    # Compute chromosome size
    group_by(CHR) %>% 
    summarise(chr_len=as.numeric(max(BP))) %>% 
    
    # Calculate cumulative position of each chromosome
    mutate(tot=cumsum(chr_len)-chr_len) %>%
    select(-chr_len) %>%
    
    # Add this info to the initial dataset
    left_join(gwasResults, ., by=c("CHR"="CHR")) %>%
    
    # Add a cumulative position of each SNP
    arrange(CHR, BP) %>%
    mutate( BPcum=BP+tot) %>%
    
    # Add highlight and annotation information
    mutate( is_highlight=ifelse(snp_trait %in% snp_trait_toAnnotate, "yes", "no")) %>%
    mutate( is_annotate=ifelse(snp_trait %in% snp_trait_toAnnotate, "yes", "no")) 
  
  print('Prepare X axis')
  # Prepare X axis
  axisdf <- don %>% group_by(CHR) %>% summarize(center=( max(BPcum) + min(BPcum) ) / 2 )
  
  print('Make the plot')
  # Make the plot
  tol21rainbow= c("#771155","#114477","#117777","#117744","#777711","#774411","#771122", 
                  "#AA4488","#4477AA", "#44AAAA","#44AA77", "#AAAA44","#AA7744","#AA4455",
                  "#CC99BB", "#77AADD", "#77CCCC", "#88CCAA", "#DDDD77", "#DDAA77", "#DD7788",
                  "#117733")
  tol7rainbow = rep(c("#CC99BB", "#77AADD", "#77CCCC", "#88CCAA", "#DDDD77", "#DDAA77", "#DD7788"),4)[1:22]
  
  ####Start subsampling (only needed for saving as PDF)
  ####Here comes a part where I reduces the filesize of the pdf by subsampling the many overlapping points
  #For the very unsignificant points, subsample so as to not have so many points
  #Divide the section between the max p-value and the first cutoff into X sections and reduce the numbe of points
  #in each section to a certain maximum. 
  #Sort by pval
  don$indices = order(don$P, decreasing = T)
  don$minLog10P = -log10(don$P)
  #Only remove points that are below the least significant pval cutoff
  toReduceDon = which(don$P>pValCutoffMark1)
  donSubset = don[toReduceDon,]
  
  #Now divide the pvalue list into X parts and 
  #subsample per section. Subsample stronger
  #From higher p-values, since these contain more SNPs
  #Always keep a minimum number of SNPs in each section, so there
  #Will still be enough points to fill the area//
  startPval = round(max(don$P), digits=6)
  startPvalLog = -log10(startPval)
  pValCutoffMark1Log = -log10(pValCutoffMark1)
  diffMinMax = pValCutoffMark1Log - startPvalLog
  nrSteps = 200
  maxNrPerPart = nrow(donSubset)/(nrSteps*dividingFactorNonSignifSnps)
  stepSize = diffMinMax/nrSteps
  
  totalIndicesToKeep=c()
  for (stepIndex in c(1:nrSteps)){
    cutOffStart = startPvalLog + stepSize * (stepIndex-1)
    cutOffStartPlus1 = cutOffStart + stepSize
    #Check how many values are in this section
    selIndices = which(don$minLog10P<cutOffStartPlus1 & don$minLog10P>=cutOffStart)
    nrIndices = length(selIndices)
    donSubset = don[selIndices,,drop=F]
    if (nrIndices>maxNrPerPart){
      #Take out this subset
      #Select X equidistant points
      #https://stackoverflow.com/questions/23145595/sample-equidistant-points-from-a-numeric-vector
      #Keep every Nth SNP, that way the subsampling is even more eqyal
      overAllIndices = c(1:length(selIndices))
      ideal <- seq(min(overAllIndices),max(overAllIndices),(max(overAllIndices)-min(overAllIndices))/(maxNrPerPart-1))
      result <- sapply(ideal, function(x) overAllIndices[which.min(abs(overAllIndices-x))] )
      toKeepIndices = selIndices[result]
    }else{
      toKeepIndices = selIndices
    }
    totalIndicesToKeep = c(totalIndicesToKeep,toKeepIndices)
    
  }
  
  signifDon = which(don$P<=pValCutoffMark1)
  toKeepTotal =   c(totalIndicesToKeep,signifDon)
  don = don[toKeepTotal,]
  #####END OF SUBSAMPLING
  
  #Define some colors
  tolGreys = rep(c( "#8c8c8c", "#b7b7b7"),11)
  names(tolGreys) = c(1:22)
  
  #Color the cytokines by stimulus and color the circulating markers with shades of blue
  if (colorByStims){
    don$stimShort = as.factor(mapLongCytoToShortCytoName(str_split_fixed(don$trait, pattern = '_', 2)[,1]))
    colorsSelStims = c('#55a868',"#dd8452","#8172b3","#da8bc3")
    names(colorsSelStims) = c("Mt","Sa","Ca","LPS")
    colorSet = c(tolGreys,colorsSelStims[levels(don$stimShort)])
    don$trait = gsub('\\.fc$','',don$trait)
    don$cyto = str_split_fixed(don$trait, pattern = '_', 4)[,4]
    don$cytoMapped = unlist(lapply( don$cyto, function(x) mappingListCyto[[x]]))
    don$snp_trait_short = paste(don$cytoMapped,paste0('(',don$stimShort,')'),don$SNP)
  }else{
    don$stimShort = as.factor(paste0('chr',don$CHR))
    colorsCircMark = rep(c( "#4c72b0", "#769ad1"),11)
    names(colorsCircMark) = levels(don$stimShort)
    colorSet=c(tolGreys,colorsCircMark)
    don$cyto = str_split_fixed(don$trait, pattern = '_', 2)[,2]
    don$cytoMapped = unlist(lapply( don$cyto, function(x) mappingListCircMed[[x]]))
    don$snp_trait_short = paste(don$cytoMapped,don$SNP)
  }
  
  don_nonSignif = don[which(don$P>pValCutoffMark1),]
  don_signif = don[which(don$P<=pValCutoffMark1),]
  
  # custom X axis:
  customLab = c("1","2","3","4","5","6","7","8","9","10","","12","","14","","16","","18","","","","22")
  
  #Make the plot
  p = ggplot(don,  aes(x=BPcum, y=-log10(P))) +
    # Show all points
    #First plot the non-significant ones
    changePlotRasterization(geom_point(data = don_nonSignif,  aes(x=BPcum, y=-log10(P), color=as.factor(CHR)), alpha=0.8, size=1.3), rasterizePlot=rasterizePlot) +
    #Plot the significant ones
    changePlotRasterization(geom_point(data = don_signif,  aes(x=BPcum, y=-log10(P), color=stimShort), alpha=0.8, size=1.3), rasterizePlot=rasterizePlot) +
    scale_color_manual(values = colorSet) +
    
    
    #scale_x_continuous( label = axisdf$CHR, breaks= axisdf$center, expand = expansion(mult = c(0.008,0.008)) ) +
    scale_x_continuous( label = customLab, breaks= axisdf$center, expand = expansion(mult = c(0.008,0.008)) ) +
    scale_y_continuous(expand = expansion(mult = c(0,0.05) )) +     # remove space between plot area and x axis
    
    if(highlightSignificant){  
      # Add highlighted points
      p = p + changePlotRasterization(geom_point(data=subset(don, is_highlight=="yes"), color="#ccb974", size=2) , rasterizePlot=rasterizePlot)
    }
  p = p + geom_hline(yintercept = -log10(pValCutoffMark1), linetype= "dashed", color = "#4c72b0") 
  p = p + geom_hline(yintercept = -log10(pValCutoffMark2), linetype= "dashed", color = "#4c72b0") 
  # Add label using ggrepel to avoid overlapping
  p = p + 
    #To create a partial transparancy, plot double
    geom_label_repel( data=subset(don, is_annotate =="yes"), aes(label=snp_trait_short), ylim  = c(6.5,NA), nudge_y=1.5, 
                      seed=1234, alpha = 0.5, force=2, segment.alpha=.5, label.padding=.1, label.size = NA, size = textSize * (5/14.2256)) +
    geom_label_repel( data=subset(don, is_annotate=="yes"), aes(label=snp_trait_short), ylim  = c(6.5,NA), nudge_y=1.5, 
                      seed=1234, fill = NA, alpha = 1, force=2, segment.alpha=.5, label.padding=.1, label.size = NA, size = textSize * (5/14.2256)) +
    ylab(expression(-log[10]*italic(P)*textstyle("-")*value)) + 
    # Custom the theme:
    theme_bw() +
    theme( 
      legend.position="none",
      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank(),
      text = element_text(family = "Helvetica", colour = 'black'),
      axis.text = element_text(size=textSize, color='black'),
      axis.title = element_text(size=textSize, color='black')) +
    xlab('Chromosome')
  
  return(p)
}


convertQtlToCorrectFormatForPlot <- function(qtlDf){
  print("select columns")
  gwasResults <- data.frame(SNP= qtlDf$snps,  
                            CHR= as.integer(gsub('chr','',qtlDf$CHROM)),
                            BP= qtlDf$POS,
                            P= qtlDf$pvalue,
                            trait = qtlDf$cytokine)
  
  #I do not want to use a hard p-value cut-off any more
  #First prune the SNPs down until we have 10 that are at least 1MB apart
  gwasResultsOfInterest = gwasResults
  gwasResultsOfInterest = gwasResultsOfInterest[order(gwasResultsOfInterest$P),,drop=F]
  return(gwasResultsOfInterest)
}


defineSnpsToPotentiallyAnnotate <- function(gwasResultsOfInterest,minNrOfSnpsToAnnotate, pValCutoffToAnotate, distSel=1e6){    
  #Based on just 1MB distances determine the p value to annotate 10 independent loci
  chromPosList = lapply(c(1:22), function(x) c())
  names(chromPosList) = c(1:22)
  snpIndex = 0
  while (length(unlist(chromPosList))<minNrOfSnpsToAnnotate) {
    snpIndex = snpIndex + 1
    print(snpIndex)
    chrSel = gwasResultsOfInterest[[snpIndex,'CHR']]
    posSel = gwasResultsOfInterest[[snpIndex,'BP']]
    minDist = min(abs(chromPosList[[chrSel]] -  posSel),na.rm=T)
    print(minDist)
    #only annotate another SNP if it is al least X distance away
    if (length(chromPosList[[chrSel]])==0 || (minDist  > distSel  )){
      chromPosList[[chrSel]] = c(chromPosList[[chrSel]], posSel)
    }
    currentPVal = gwasResultsOfInterest[[snpIndex,'P']]
    print(currentPVal)
  }
  
  #Take the maximum p-alue from either the top 10 SNPs or the p-value cutoff defined manually
  if (currentPVal>pValCutoffToAnotate){
    pValCutoffToAnotateFinal = currentPVal
  }else{
    pValCutoffToAnotateFinal = pValCutoffToAnotate
  }
  
  gwasResultsOfInterest = gwasResultsOfInterest[gwasResultsOfInterest$P<= pValCutoffToAnotateFinal,,drop=F]
  
  #Iterate over the SNPs, and for the SNPs in linkage take the top snp
  #Otherwise the plot is too crowded
  print('Iterate over the SNPs, and take just SNPs further than 1MB on each side or from different quantitative trait, if in linkage')
  gwasResultsOfInterestFilt = gwasResultsOfInterest
  print(head(gwasResultsOfInterestFilt))
  snpIndex = 0
  nrOfSnp = 0
  
  #Further remove SNPs from the same cyt-stimulus within a certain area
  toRemoveIndices = c()
  for (snpIndex in 1:nrow(gwasResultsOfInterest)){
    print(paste0(snpIndex,' out of ',nrow(gwasResultsOfInterest)))
    if (snpIndex %in% toRemoveIndices){next}
    
    print(snpIndex)
    traitSel = as.character(gwasResultsOfInterest[[snpIndex,'trait']])
    print(traitSel)
    rsIdSel = as.character(gwasResultsOfInterest[[snpIndex,'SNP']])
    print(rsIdSel)
    
    rowsForThisStimulus = which(gwasResultsOfInterest$trait==traitSel)
    gwasResultsOfInterestTemp = gwasResultsOfInterest[rowsForThisStimulus,,drop=F]
    
    #Remove all snps for this stimulus within 100kB 
    distanceWithSnps = gwasResultsOfInterest[[snpIndex,'BP']] - gwasResultsOfInterestTemp[,'BP']
    
    toRemove = gwasResultsOfInterestTemp[which(abs(distanceWithSnps)<distSel),'SNP']
    toRemove = setdiff(toRemove,rsIdSel)
    toRemoveRows = (gwasResultsOfInterest$SNP %in% toRemove) & gwasResultsOfInterest$trait==traitSel
    toKeepRows = !toRemoveRows
    toRemoveIndices = c(toRemoveIndices,which(toRemoveRows))
    
  }
  gwasResultsOfInterestFilt = gwasResultsOfInterest[setdiff(c(1:nrow(gwasResultsOfInterest)),toRemoveIndices),,drop=F]
  #Return a filtered df with just the ones we  want to annotates
  return(gwasResultsOfInterestFilt)
}


source('processingFxnsGeneral.R')
source('plottingFxnsGeneral.R')
source('300bcg_generalFxns.R')

####EXAMPLES OF HOW TO RUN#####
# Rscript 300bcg_topCQtlHits.R cytoData.good.log2 V1
# Rscript 300bcg_topCQtlHits.R cytoData.good.fc V3
# Rscript 300bcg_topCQtlHits.R circMarkerData.log2 V1
# Rscript 300bcg_topCQtlHits.R circMarkerData.fc V3
# sbatch 300bcg_topCQtlHits.batch cytoData.good.fc V3
# sbatch 300bcg_topCQtlHits.batch cytoData.good.fc V2
# sbatch 300bcg_topCQtlHits.batch circMarkerData.log2 V1
# sbatch 300bcg_topCQtlHits.batch cytoData.good.log2 V1
#args=c("cytoData.good.log2", "V1")
#args=c("cytoData.good.fc", "V2")
#args=c("cytoData.good.fc", "V3")
#args=c("circMarkerData.log2", "V1")

##Load some mapping lists that define how to translate the short names to long names
mapCytokineNamesFile = file.path('StimulusToPrettyNames.xlsx')
sheetName = 'cytokine'
loadedData = openxlsx::read.xlsx(xlsxFile = mapCytokineNamesFile, sheet = sheetName, colNames=F)
mappingListCyto = createMappingList(loadedData)

sheetName = 'circMed'
loadedData = openxlsx::read.xlsx(xlsxFile = mapCytokineNamesFile, sheet = sheetName, colNames=F)
mappingListCircMed = createMappingList(loadedData)

#Get the system arguments
args = base::commandArgs(trailingOnly=TRUE)
print(paste0('arguments=',args))
print(length(args))

#Load cytokine and metadata
allRelevantDataList = loadLukasMetaDataPbmc(allMetaDataPath)
list2env(allRelevantDataList, .GlobalEnv)

################SOME SETTINGS################
#Define for which data to run the algorithm 
dataNameSel = args[[1]]
nameToAddToFolder = dataNameSel
selData = eval(parse(text = dataNameSel))
selVariables = colnames(selData) 

mainTpSel = args[[2]]
splitSetting="perChrom"

###Get the variables we need
library(config)
Sys.setenv(R_CONFIG_ACTIVE = "300bcg")
config <- config::get()
for (x in names(config)){eval(parse(text=paste0(x,"='",config[[x]],"'")))}

####Set a value that is used in the Manhattan plot
#The problem with the Manhattan is that is will become quite large (file size)
#due to the large number of points plotted
#Many of them will overlap though, so some can be removed
#In the manhattan function I included some code that gets rid of points in a 
#somewhat intelligent way. However, this needs some parameters to be tuned
#For now I hardcoded it below here
#The value below should be roughly the maximum -log10(pvalue) in the plot
if (dataNameSel=="cytoData.good.log2" & mainTpSel=='V1'){
  dividingFactorNonSignifSnps = 15
}else if (dataNameSel=="circMarkerData.log2" & mainTpSel=='V1'){
  dividingFactorNonSignifSnps = 37
}else if (dataNameSel=="cytoData.good.fc" & (mainTpSel %in% c('V3','V2'))){
  dividingFactorNonSignifSnps = 7
}else if (dataNameSel=="circMarkerData.fc" & (mainTpSel%in% c('V3','V2'))){
  dividingFactorNonSignifSnps = 5
}

SuppTableExtendedMapping = list()
SuppTableExtendedMapping[["circMarkerData.log2_V1"]] = "Ext_S05.1"
SuppTableExtendedMapping[["cytoData.good.log2_V1"]] = "Ext_S05.2"
SuppTableExtendedMapping[["cytoData.good.fc_V2"]] = "Ext_S05.3"
SuppTableExtendedMapping[["cytoData.good.fc_V3"]] = "Ext_S05.4"

#####Correct the data for top hits plotting####
#We also want to plot some examples. 
#Define the covariates
correctionFactorList = loadCorrectionFactorList()
covariatesSel = correctionFactorList[[dataNameSel]]

#Define individuals to remove
if (endsWith(dataNameSel, '.fc')){
  #For the fold-changes we remove the people vaccinated in the evening
  toRemoveEve = rownames(allRelevantData)[allRelevantData$IC_EVENING]
}else{
  toRemoveEve = NULL
}

#Outliers
#Define who to exclude based on some other parameters
outliers = rownames(metaData)[which(rowSums(metaData[,c("EXCLUSION","SNP_OUTLIER")]!='')>0)]

#Keep just those indiv combined with the evening in those cases where we wanted to remove them for other reasons
toRemoveIndividuals = c(outliers,toRemoveEve)

#Remove 
allRelevantData = allRelevantData[which(!(rownames(allRelevantData) %in% toRemoveIndividuals)),]
subFolderForSave = paste('quantitative',dataNameSel,paste(mainTpSel,collapse='.'),sep='__')

###Define some pretty name to use when saving 
prettyNameForFileName = paste0('SNP_',definitionsForSupplement[[dataNameSel]],'_',definitionsForSupplement[[mainTpSel]])

#check which results are there
qtlDirOverall = file.path(genetics_save_dir,"QTLresults")
qtlDir = file.path(qtlDirOverall,subFolderForSave,splitSetting)
allQtlFiles = list.files(qtlDir, pattern = ".*_withInfo.cqtls")
print('allQtlFiles')
print(allQtlFiles)
print('selVariables')
print(selVariables)
allCytosInQtlFiles = unique(gsub(pattern = '(\\_chr[0-9]{1,2})?\\_withInfo\\.cqtls$',replacement = '',allQtlFiles))
print('allCytosInQtlFiles')
print(allCytosInQtlFiles)

##Check for which cytokines we also have results
allCytosInQtlFiles = intersect(selVariables,allCytosInQtlFiles)

dirToSaveResultsQtlCombine = file.path(dirname(qtlDir),'combinedResults')
mkdirs(dirToSaveResultsQtlCombine)

###We also want a supplemental table showing all hits, lets combine all this data in a smart way
###I do not want to leave out any SNPs in the online FAIR database
saveSupplementalTable = F
if (saveSupplementalTable){
  
  ##Load variant info
  print('combining variantinfo')
  variantInfoCombined = NULL
  for (chrSel in c(1:22)){
    print(chrSel)
    subFolderWithData='mafHweOutlierFilteredGrCh38'
    rdsOutput=file.path(folderWithGenotypeData,'RDS',subFolderWithData)
    
    saveFileGeneticsChrom = file.path(rdsOutput,paste0("GRCh38_300BCG_chr_",chrSel,"_annotated_MAF0.1_HWE1e-5_noOutliers.recode.tags_dosages.RDS"))
    variantInfo = readRDS(file.path(rdsOutput,paste0("chr_",chrSel,"_300bcg_variantInfoAddedChrPos.RDS")))
    variantInfo.sel = variantInfo[,c("ID","CHROM","POS","REF","ALT")]
    colnames(variantInfo.sel) = c(definitionsForSupplement[['snp']],definitionsForSupplement[['chr']],definitionsForSupplement[['pos']],
                                  definitionsForSupplement[['ref']],definitionsForSupplement[['alt']])
    if (is.null(variantInfoCombined)){
      variantInfoCombined = variantInfo.sel
    }else{
      variantInfoCombined = rbind.data.frame(variantInfoCombined,variantInfo.sel)
    }
  }
  print("variantInfoCombined")
  print(head(variantInfoCombined))
  
  sortedCols = c(definitionsForSupplement[["snp"]],definitionsForSupplement[["chr"]],definitionsForSupplement[["pos"]],
                 definitionsForSupplement[["ref"]],definitionsForSupplement[["alt"]],
                 definitionsForSupplement[[paste(dataNameSel,'LONG',sep='_')]],definitionsForSupplement[["coef"]],
                 definitionsForSupplement[["tval"]],definitionsForSupplement[["pval"]])
  
  partialTable = NULL
  for (cytokineSel in setdiff(colnames(allRelevantDataList[[dataNameSel]]),cytosToExclude)){
    ##First go over all chromosomes and create a full set per cytokine
    fullPerCytoDf = NULL
    for (chrSel in c(1:22)){
      fileSel = paste0(cytokineSel,"_chr",chrSel,"_withInfo.cqtls")
      
      print(fileSel)
      qtlRes = fread(file.path(qtlDir,fileSel))
      qtlRes = qtlRes[,c("snps","gene","beta","statistic","pvalue")]
      colnames(qtlRes) = c(definitionsForSupplement[['snp']],
                           definitionsForSupplement[[paste(dataNameSel,'LONG',sep='_')]],
                           definitionsForSupplement[['coefficient']],
                           definitionsForSupplement[['tval']],
                           definitionsForSupplement[['pval']])
      selColSnps = definitionsForSupplement[['snp']]
      qtlRes = qtlRes[!is.na(unlist(qtlRes[,..selColSnps])),]
      #Remove any for which the SNPs were not assigned a name or chrom pos location
      qtlRes = qtlRes[which(unlist(qtlRes[,..selColSnps])!=""),]
      
      if (is.null(fullPerCytoDf)){
        fullPerCytoDf = qtlRes
      }else{
        fullPerCytoDf = rbind.data.frame(fullPerCytoDf,qtlRes)
      }
    }
    
    ###Now map the cytokine to prettier name
    factorCol = definitionsForSupplement[[paste(dataNameSel,'LONG',sep='_')]]
    print('substitute')
    cytoPart = gsub('\\_good$','',cytokineSel)
    cytoPart = gsub('\\.fc$','',cytokineSel)
    print(cytoPart)
    print('map')
    cytoPart.map = mapCytoNamesForSupplement(cytoPart,dataNameSel)
    print(cytoPart.map)
    print('replace in df')
    fullPerCytoDf[,factorCol] = cytoPart.map
    
    print("fullPerCytoDf")
    print(head(fullPerCytoDf))
    
    # ####add variant info to the table
    print('merge results with variantinfo')
    snpCol = definitionsForSupplement[['snp']]
    fullPerCytoDf = base::merge(fullPerCytoDf,variantInfoCombined, by=snpCol, all.x=T, all.y=F)
    #Sort
    pValCol = definitionsForSupplement[['pval']]
    fullPerCytoDf = fullPerCytoDf[order(unlist(fullPerCytoDf[,..pValCol]), decreasing = F),]
    
    #Now save the results as csv
    folderToSaveCsv = file.path(suppTableDir,'csv_large',paste0(SuppTableExtendedMapping[[paste(dataNameSel,mainTpSel,sep='_')]],"_GWAS_",prettyNameForFileName))
    
    #Clean names and round
    cleanCytoNameForFilename=janitor::make_clean_names(cytoPart.map,case='none')
    mkdirs(folderToSaveCsv)
    fullPerCytoDf[[definitionsForSupplement[['tval']]]] = round(fullPerCytoDf[[definitionsForSupplement[['tval']]]], digits = 4)
    fullPerCytoDf[[definitionsForSupplement[['coef']]]] = round(fullPerCytoDf[[definitionsForSupplement[['coef']]]], digits = 4)
    fullPerCytoDf[[definitionsForSupplement[['pval']]]] = signif(fullPerCytoDf[[definitionsForSupplement[['pval']]]], digits = 4)
    
    fullPerCytoDf = fullPerCytoDf[,..sortedCols]
    
    fwrite(x = fullPerCytoDf, quote = F, file = file.path(folderToSaveCsv,paste0(cleanCytoNameForFilename,'_',prettyNameForFileName,'_GWAS.csv')), row.names = F, col.names = T, sep=';')
    
    #Take the significant ones with to a max of X lines, 'fullTable" has to be sorted by pvalue
    toSave = intersect(c(1:5000),which(fullPerCytoDf[,..pValCol]<1e-4))
    partialTablePerCyto = fullPerCytoDf[toSave,]
    
    ##Now merge all data for the main
    if (is.null(partialTable)){
      partialTable = partialTablePerCyto
    }else{
      partialTable = rbind.data.frame(partialTable,partialTablePerCyto)
    }
  }
  
  
  print('dimensions of partialTable')
  print(dim(partialTable))
  #Sort rows
  partialTable = partialTable[order(unlist(partialTable[,..pValCol]), decreasing = F),]
  
  #Sort columns
  partialTable = partialTable[,..sortedCols]
  
  #For the actual supplement we are going to use xlsx in the end
  #We can use multiple tabs for this
  #Use "subFolderForSave" as a tab
  #Add to a file if it already exists
  #From https://github.com/awalker89/openxlsx/issues/157
  
  folderToSaveXlsx = file.path(suppTableDir,'xlsx_small','raw')
  mkdirs(folderToSaveXlsx)
  fileNameXlsx = file.path(folderToSaveXlsx,paste0('TableS5_',tableAdditionToFileNames[["TableS5"]],'.xlsx'))

  if (!(file.exists(fileNameXlsx))){
    wb <- openxlsx::createWorkbook()
  }else{
    wb <- openxlsx::loadWorkbook(fileNameXlsx)
  }
  
  sheetName = prettyNameForFileName
  
  if (!(sheetName %in% wb$sheet_names)){
    openxlsx::addWorksheet(wb = wb, sheet = sheetName)
  }
  
  bold_style <- createStyle(textDecoration = "Bold")
  # make non-numeric column bold
  selCols = which(unlist(lapply(partialTable,class))!="numeric")
  
  openxlsx::addStyle(wb, style = bold_style, rows = 1:(nrow(partialTable)+1), sheet = sheetName, cols = selCols, gridExpand = T)
  
  # create style, in this case bold header
  header_st <- openxlsx::createStyle(textDecoration = "Bold")
  
  openxlsx::writeData(wb = wb, sheet = sheetName, x = partialTable, colNames = T, rowNames = F, headerStyle = header_st)
  openxlsx::saveWorkbook(wb = wb, file = fileNameXlsx, overwrite = T)
  
}



createPlots = T
if (createPlots){
  totalQtlRes = loadTopHitsQtl(dirToSaveResultsQtlCombine,qtlDir, allQtlFiles, cytosToExclude = NULL, reExtractTopHits=F )
  
  #Add chrom pos for those who have it missing --> this was a problem with how I originally ran the analysis
  #There were samples with missing RS-ids. When adding the chrom as pos back into the results files
  #it could not find a match with the original names for those without rs-ids. Luckily this info is contained in
  #the SNP names for these cases
  noChromPos = which(startsWith(totalQtlRes$snps, prefix = 'chr'))
  snpsToChange = unlist(totalQtlRes[noChromPos,c('snps')])
  if (length(snpsToChange)>0){
    snpsToChange.split = t(as.data.frame(strsplit(snpsToChange, "\\:|\\_|\\/")))
    totalQtlRes[noChromPos,'CHROM'] = snpsToChange.split[,1]
    totalQtlRes[noChromPos,'POS'] = as.numeric(snpsToChange.split[,2])
    totalQtlRes[noChromPos,'REF'] = snpsToChange.split[,3]
    totalQtlRes[noChromPos,'ALT'] = snpsToChange.split[,4]
  }
  
  ####Create manhattan plot for the top hits
  #Define two cut-off threshold to show
  pValCutoffMark1 = 5e-6
  pValCutoffMark2= 5e-8
  #Convert to a format the plot function likes
  gwasResultsOfInterest = convertQtlToCorrectFormatForPlot(totalQtlRes)
  #Define which SNPs to annotate potentially
  minNrOfSnpsToAnnotate=5
  pValCutoffToAnotate=5e-8
  snpsToAnnotate = defineSnpsToPotentiallyAnnotate(gwasResultsOfInterest,minNrOfSnpsToAnnotate,pValCutoffToAnotate)
  #Create manhattan plot]
  
  if (startsWith(x = dataNameSel, prefix = 'cytoData')){
    colorByStims = T
  }else{
    colorByStims = F
  }
  
  dirToSaveQtlImages = file.path(image_dir_300BCG,'manuscriptImagesRaw','genetics','manhattan',prettyNameForFileName)
  mkdirs(dirToSaveQtlImages)
  
  colorOutline ='#937860'
  p = createManhattanPlot(gwasResultsOfInterest, snpsToAnnotate, pValCutoffMark1, pValCutoffMark2, colorByStims = colorByStims,textSize=12,
                          dividingFactorNonSignifSnps= dividingFactorNonSignifSnps,colorOutline = colorOutline,rasterizePlot = T)
  
  if (dataNameSel=="circMarkerData.log2"){
    widthSel=5
    widthPanel = 2*5.1216*0.393701
    heightSel = 5
    heightPanel = 4.47
  }else{
    #Set the panel size for this
    widthSel=5
    widthPanel = 2*5.1216*0.393701#4.47
    heightSel = 5
    heightPanel = 4.47
  }
  
  library(egg)
  #Fix the dimensions of the panel, so all of them nicely align
  p_fixed <- set_panel_size(p,
                            width  = unit(widthPanel, "in"),
                            height = unit(heightPanel, "in"))
  
  fileToSaveCairo = file.path(dirToSaveQtlImages,paste0('manhattanPlot','_',dataNameSel,'_',mainTpSel,'_cairo.pdf'))
  ggsave(filename = fileToSaveCairo, plot = p_fixed, width = widthSel, height=heightSel, device = cairo_pdf)
  
  fileToSaveCairo2 = file.path(dirToSaveQtlImages,paste0('manhattanPlot','_',dataNameSel,'_',mainTpSel,'_cairo2.pdf'))
  cairo_pdf(fileToSaveCairo2, family = "Helvetica")
  print(p_fixed)
  dev.off()
  
  fileToSave = file.path(dirToSaveQtlImages,paste0('manhattanPlot','_',dataNameSel,'_',mainTpSel,'.pdf'))
  ggsave(filename = fileToSave, plot = p_fixed, width = widthSel, height=heightSel, useDingbats=F)
  print(paste('saved file',fileToSave))
  
  ######Also plot the top hits themselves
  plotTopHits=T
  examplePlotsFolder = file.path(image_dir_300BCG,'manuscriptImagesRaw','genetics','examplePlots',prettyNameForFileName)
  mkdirs(examplePlotsFolder)
  #Visualising the top hits from the genetics
  if (plotTopHits){
    totalQtlRes = as.data.frame(totalQtlRes)
    totalQtlResToPlot = totalQtlRes[which(totalQtlRes$snps %in% snpsToAnnotate$SNP),]
    for (rowIndex in c(1:nrow(totalQtlResToPlot))){
      #Check which chromosome
      chrSel = totalQtlResToPlot[rowIndex,'CHROM']
      rsIdSel = totalQtlResToPlot[rowIndex,'snps']
      cytokineSel = totalQtlResToPlot[rowIndex,'cytokine']
      #Load the relevant genetics data
      subFolderWithData='mafHweOutlierFilteredGrCh38'
      rdsOutput=file.path(folderWithGenotypeData,'RDS',subFolderWithData)
      saveFileGeneticsChrom = file.path(rdsOutput,paste0("GRCh38_300BCG_chr_",gsub('chr','',chrSel),"_annotated_MAF0.1_HWE1e-5_noOutliers.recode.tags_hardCalls.RDS"))
      
      ##Load the genetics, and copy it multiple times to each timepoint
      geneticsData = readRDS(saveFileGeneticsChrom)
      colnames(geneticsData) = paste0(colnames(geneticsData),'_V1')
      #Tp1
      geneticsDataSel = t(geneticsData[rsIdSel,,drop=F])
      geneticsDataSelTp1 = geneticsDataSel
      #Tp2
      geneticsDataSelTp2 = geneticsDataSel
      rownames(geneticsDataSelTp2) = gsub('\\_V1','\\_V2',rownames(geneticsDataSelTp2))
      #Tp3
      geneticsDataSelTp3 = geneticsDataSel
      rownames(geneticsDataSelTp3) = gsub('\\_V1','\\_V3',rownames(geneticsDataSelTp3))
      
      geneticsDataSel = rbind(geneticsDataSelTp1,geneticsDataSelTp2,geneticsDataSelTp3)
      
      #Map to reference and alternative
      #Load variantinfo
      fileWithAnno = file.path(rdsOutput,paste0("chr_",gsub('chr','',chrSel),"_300bcg_variantInfoAddedChrPos.RDS"))
      annotationInfo = readRDS(fileWithAnno)
      rowIndexAnnotation = which(annotationInfo$POS==totalQtlResToPlot[rowIndex,'POS'])
      refAllele = as.character(annotationInfo[[rowIndexAnnotation,'REF']])
      altAllele = as.character(annotationInfo[[rowIndexAnnotation,'ALT']])
      
      #Convert numeric variant to character
      geneticsDataSel[,rsIdSel] = as.character(geneticsDataSel[,rsIdSel])
      geneticsDataSel[which(geneticsDataSel[,rsIdSel]=="0"),] = paste0(refAllele,refAllele)
      geneticsDataSel[which(geneticsDataSel[,rsIdSel]=="1"),] = paste0(refAllele,altAllele)
      geneticsDataSel[which(geneticsDataSel[,rsIdSel]=="2"),] = paste0(altAllele,altAllele)
      
      ##Now for the plot regress out the covariates
      #Linear model correction
      combinedDataSelAll = base::merge(geneticsDataSel,allRelevantData, by='row.names')
      rownames(combinedDataSelAll) = combinedDataSelAll$Row.names
      combinedDataSelAll$Row.names <- NULL
      #Remove data which is not missing in the selected data
      combinedDataSelAll = combinedDataSelAll[which(!is.na(combinedDataSelAll[,cytokineSel])),]
      
      combinedDataSelAll[, rsIdSel] <- factor(combinedDataSelAll[, rsIdSel], levels = c(paste0(refAllele,refAllele),paste0(refAllele,altAllele),paste0(altAllele,altAllele)))
      
      if (endsWith(x = dataNameSel, suffix = 'fc')){
        timepointsToPlot = c('V2','V3')
      }else{
        timepointsToPlot = c('V1','V2','V3')
      }
      
      combinedDataSelAll.cor = combinedDataSelAll
      for (timepointSel in timepointsToPlot){
        indicesWithTimepoint = endsWith(timepointSel, x = rownames(combinedDataSelAll))
        dataSubsetSel = combinedDataSelAll[indicesWithTimepoint,]
        formulaSel = paste0(cytokineSel,' ~ ',paste(covariatesSel, collapse = ' + '))
        linearFit = lm(formula = formulaSel, data = dataSubsetSel)
        residuals.sel = resid(linearFit) + mean(dataSubsetSel[,cytokineSel],na.rm=T)
        combinedDataSelAll.cor[names(residuals.sel),cytokineSel] = residuals.sel
      }
      
      library(reshape2)
      cytokineSel.full = gsub('\\.fc$','',cytokineSel)
      if (colorByStims){
        shortStimSel = mapLongCytoToShortCytoName(str_split_fixed(cytokineSel.full, pattern = '_', 2)[,1])
        cytoMapped = mappingListCyto[[str_split_fixed(cytokineSel.full, pattern = '_', 4)[,4]]]
        cytoTotal = paste0(cytoMapped, ' (', shortStimSel,')')
        titleSel = paste(cytoTotal,rsIdSel)
        
        ylabSel = cytoTotal
      }else{
        cytoMapped = mappingListCircMed[[str_split_fixed(cytokineSel.full, pattern = '_', 2)[,2]]]
        titleSel = paste(cytoMapped,rsIdSel)
        cytoTotal = cytoMapped
        ylabSel = cytoTotal
      }
      
      #Now make the plots
      toPlotSingle = combinedDataSelAll.cor[which(endsWith(suffix = mainTpSel,x = rownames(combinedDataSelAll.cor))),]
      
      #Now transform back if needed
      if (dataNameSel=="cytoData.good.log2"){
        toPlotSingle[,cytokineSel] = 2^toPlotSingle[,cytokineSel]
      }
      
      set.seed(10)
      p = ggplot(toPlotSingle, aes_string(rsIdSel, cytokineSel, fill=rsIdSel)) +
        geom_boxplot(width=0.7, outlier.shape = NA) + geom_jitter(width = .3, height=0) +
        theme_bw() +
        theme( 
          legend.position="none",
          panel.grid.minor = element_blank(),
          text = element_text(family = "Helvetica", colour = 'black'),
          axis.text = element_text(size=12, color='black'),
          axis.title = element_text(size=12, color='black'),
          plot.title = element_text(hjust = 0.5, color='black')
        ) + 
        scale_fill_manual(values = c("#4c72b1", "#769ad1","#9bbbea")) +
        ylab(ylabSel) + xlab(rsIdSel)
      
      if (dataNameSel=="cytoData.good.log2"){
        p = p + scale_y_continuous(trans = scales::log_trans(base = 10), breaks = scales::pretty_breaks(n=10))#scales::log_breaks(base=10))
      }
      
      print(p)
      
      sizeCombis = list()
      sizeCombis[["1"]] = c(2, 1.2, 5, 4.47)
      sizeCombis[["2"]] = c(2, 1.2, 2.3, 2)
      
      for (indexSize in names(sizeCombis)){
        
        widthSel= sizeCombis[[indexSize]][[1]]
        widthPanel = sizeCombis[[indexSize]][[2]]
        heightSel = sizeCombis[[indexSize]][[3]]
        heightPanel = sizeCombis[[indexSize]][[4]]
        
        library(egg)
        #Fix the dimensions of the panel, so all of them nicely align
        p_fixed <- set_panel_size(p,
                                  width  = unit(widthPanel, "in"),
                                  height = unit(heightPanel, "in")
        )
        
        ggsave(p_fixed, filename = file.path(examplePlotsFolder,paste0('examplePlot_',rsIdSel,'_',cytokineSel,'_covariateCorrectedSingleTp_',indexSize,'.pdf')), 
               width = widthSel, height=heightSel, device = cairo_pdf)
      }
    }
  }
}