library(data.table)
library(dplyr)
#library(formattable)
library(tidyr)
library("R.utils")
library('openxlsx')
library(ggplot2)
library(stringr)
library(data.table)
library(dplyr)
library(RColorBrewer)

print('starting')

###Get the variables we need
library(config)
Sys.setenv(R_CONFIG_ACTIVE = "300bcg")
config <- config::get()
for (x in names(config)){eval(parse(text=paste0(x,"='",config[[x]],"'")))}

source(file.path('300bcg_generalFxns.R'))

######Load cytokine and metadata######
allRelevantDataList = loadLukasMetaDataPbmc(allMetaDataPath)
list2env(allRelevantDataList, .GlobalEnv)

#Load mapping lists to map to simpler cytokine names
mapCytokineNamesFile = file.path('StimulusToPrettyNames.xlsx')
sheetName = 'cytokine'
loadedData = openxlsx::read.xlsx(xlsxFile = mapCytokineNamesFile, sheet = sheetName, colNames=F)
mappingListCyto = createMappingList(loadedData)

sheetName = 'circMed'
loadedData = openxlsx::read.xlsx(xlsxFile = mapCytokineNamesFile, sheet = sheetName, colNames=F)
mappingListCircMed = createMappingList(loadedData)

#Load peak annotations
peakAnnotations = fread(peakInfoFile, data.table = F)
cytokineGroups = defineSubgroupsCytokines()

#Define how many of the hits to plots
nrOfHitsToPlot = 10

#For which combinations to make the plots
combinationsToCheck = t(as.data.frame(list(
  c("cytoData.good.fc","V3"),
  c("cytoData.good.fc","V2"),
  c("circMarkerData.log2","V1"),
  c("cytoData.good.log2","V1")
)))

#Iterate over there combis
print(combinationsToCheck)
for (rowIndex in c(1:nrow(combinationsToCheck))){
  
  dataTypeSel=combinationsToCheck[[rowIndex,1]]
  timepointSel=combinationsToCheck[[rowIndex,2]]
  
  dataSel = allRelevantDataList[[dataTypeSel]]
  selVariables = colnames(dataSel)
  
  if (startsWith(x = dataTypeSel, prefix = 'cytoData')){
    selData = 'cytokines'
  }else{
    selData = 'CM'
  }
  
  #Define folder for fold-changes and absolute data
  if (endsWith(x = dataTypeSel, suffix = '.fc')){
    selFolder = file.path(folderWithTrainingAtacResults,selData)
  }else{
    selFolder = file.path(folderWithBaselineAtacResults,selData)
  }
  
  combinedResults <- loadAllAtacSeqAssociationResults(selFolder)
  
  #Use "GENE_AND_DISTAL_10kb" region definitions for mapping regions to genes,
  for (subsetSel in c('GENE_AND_DISTAL_10kb')){
    
    #Load the gene set data
    fileWithRelevantPeaks = file.path(folderPathways,paste0('peaks_filtered_PBMC.',subsetSel,'.csv.gz'))
    peakFileSel = fread(fileWithRelevantPeaks, data.table = F)
    combinedResultsSel = combinedResults[which(combinedResults$V1 %in% peakFileSel$peak_id),]
    
    #Adjust the pvalues
    combinedResultsSel$adj.P.Val.total = p.adjust(combinedResultsSel$p.value, method = 'fdr')
    combinedResultsSel$logPVal = -log10(combinedResultsSel$p.value)
    combinedResultsSel = combinedResultsSel[order(combinedResultsSel$p.value, decreasing = F),]
    
    print(head(combinedResultsSel,nrOfHitsToPlot))
    
    #Iterate over the top hits
    for (topNr in c(1:nrOfHitsToPlot)){
      ###Now take the top hit and plot it
      #Get counts
      cytoSel = combinedResultsSel[topNr,'variable']
      peakSel = combinedResultsSel[topNr,'V1']
      
      print(paste0('cytoSel=',cytoSel))
      
      #Load the relevant data needed for the plot
      #This is based on the results from Lukas, which have a very specific file structure
      patternSel = paste0('(\\_|\\.)',timepointSel,'(\\_|\\.).*',cytoSel,'$')
      folderForCyto = list.files(selFolder, include.dirs = F, pattern = patternSel)
      
      #Load counts
      folderWithFiles = file.path(selFolder,folderForCyto)
      countsScaled = fread(file.path(folderWithFiles, paste0('de_counts_',folderForCyto,'.csv.gz')), data.table = F)
      rownames(countsScaled) = countsScaled[,1]
      countsScaled[,1] <- NULL
      
      countsScaled.sel = as.data.frame(t(countsScaled))
      countsScaled.sel$rn = rownames(countsScaled.sel)
      
      ###In the end it might be easier to get all the data from the "DESIGN" dataframe
      #This should contain all the covariates and means I do not need to do all the work to rename the column names
      folderWithFiles = file.path(selFolder,folderForCyto)
      designLoaded = fread(file.path(folderWithFiles, paste0('de_design_',folderForCyto,'.csv.gz')), data.table = F)
      rownames(designLoaded) = designLoaded[,1]
      designLoaded[,1] <- NULL
      designLoaded.sel = designLoaded
      designLoaded.sel$rn = rownames(designLoaded.sel)
      designLoaded.sel = designLoaded.sel[rownames(countsScaled.sel),]
      
      #####CORRECTION BASED ON THE COEFFICIENTS LUKAS
      #Load the coefficient data
      coefficientsSel = fread(file.path(folderWithFiles, paste0('de_coefficients_',folderForCyto,'.csv.gz')), data.table = F)
      rownames(coefficientsSel) = coefficientsSel[,1]
      coefficientsSel[,1] <- NULL
      
      cytoSelLong = colnames(coefficientsSel)[grep(pattern = paste0('.*',cytoSel,'.*'), x = colnames(coefficientsSel)) ]
      
      coefficientNamesToCorrect = setdiff(colnames(designLoaded),c("X__Intercept",cytoSelLong))
      
      #Do the actual calculation to regress out the covariate
      fittedData = as.matrix(coefficientsSel[peakSel,coefficientNamesToCorrect ]) %*% t(as.matrix(designLoaded.sel[,coefficientNamesToCorrect]))
      
      correctedData = countsScaled.sel[,peakSel,drop=F] - fittedData
      correctedData$rn = rownames(correctedData)
      
      #Merge data
      dataToPlot= base::merge(correctedData,designLoaded.sel, by='rn')
      
      #Now transform back if needed
      if (selData == 'cytokines' & endsWith(x = dataTypeSel, suffix = '.log2')){
        dataToPlot[,cytoSelLong] = 2^dataToPlot[,cytoSelLong]
      }
      
      ##Now map the cytokine name to a nicer name to use in the plot
      if (selData=='cytokines'){
        newStim = mapLongCytoToShortCytoName(str_split_fixed(cytoSel, pattern = '_', n = 5)[,1])
        newCyto = unlist(lapply(str_split_fixed(cytoSel, pattern = '_', n = 5)[,4], function(x) mappingListCyto[[x]]))
        cytoSel.short = paste0(newCyto,' (',newStim,')')
      }else if (selData=='CM'){
        newCyto = unlist(lapply(cytoSel, function(x) mappingListCircMed[[x]]))
        cytoSel.short = newCyto
      }
      
      #Select gene name
      geneSel = peakFileSel[[which(peakFileSel$peak_id==peakSel),'gene_name']]
      
      #Define color of the outline
      colorOutline ='#ccb974'
      
      #Make the plot and save
      folderToSavePlot = file.path(image_dir_300BCG,'manuscriptImagesRaw','atacSeq','examplePlots',paste(dataTypeSel,timepointSel,sep='_'))
      mkdirs(folderToSavePlot)
      
      sizeCombis = list()
      sizeCombis[["1"]] = c(2, 1.2, 5, 4.47)
      sizeCombis[["2"]] = c(2, 1.2, 2.3, 2)
      
      for (indexSize in names(sizeCombis)){
        
        widthSel= sizeCombis[[indexSize]][[1]]#2
        widthPanel = sizeCombis[[indexSize]][[2]]#1.2
        heightSel = sizeCombis[[indexSize]][[3]]#5
        heightPanel = sizeCombis[[indexSize]][[4]]#4.47
        
        p = ggplot(dataToPlot, aes_string(y=peakSel, x=cytoSelLong)) + 
          geom_point()+
          geom_smooth(method=lm, se=FALSE, color='#4c72b1') +
          theme_bw() +
          theme( 
            legend.position="none",
            text = element_text(family = "Helvetica",colour = 'black'),
            axis.text = element_text(size=12, color='black'),
            axis.title = element_text(size=12, color='black'),
            plot.title = element_text(hjust = 0.5, color='black')
          ) + 
          xlab(cytoSel.short) +
          ylab(paste0(geneSel, ' peak')) 
        library(egg)
        #Fix the dimensions of the panel, so all of them nicely align
        p_fixed <- set_panel_size(p,
                                  width  = unit(widthPanel, "in"),
                                  height = unit(heightPanel, "in"))
        
        
        ggsave(p_fixed, filename = file.path(folderToSavePlot,paste0('examplePlot_',peakSel, '_',cytoSel.short,'_covariateCorrectedSingleTp_',indexSize,'.pdf')), 
               width = 2, height=5, device = cairo_pdf)

        ###The same plot, but now with x and y flipped
        p = ggplot(dataToPlot, aes_string(x=peakSel, y=cytoSelLong)) + 
          geom_point()+
          geom_smooth(method=lm, se=FALSE, color='#4c72b1') +
          theme_bw() +
          theme( 
            legend.position="none",
            text = element_text(family = "Helvetica", colour = 'black'),
            axis.text = element_text(size=12, colour = 'black'),
            axis.title = element_text(size=12, colour = 'black'),
            plot.title = element_text(hjust = 0.5, colour = 'black')
          ) + 
          ylab(cytoSel.short) + 
          xlab(paste0(geneSel, ' peak')) 
        
        if (selData == 'cytokines' & endsWith(x = dataTypeSel, suffix = '.log2')){
          p = p + scale_y_continuous(trans = scales::log_trans(base = 10), breaks = scales::pretty_breaks(n=10))#scales::log_breaks(base=10))
        }
        
        library(egg)
        #Fix the dimensions of the panel, so all of them nicely align
        p_fixed <- set_panel_size(p,
                                  width  = unit(widthPanel, "in"),
                                  height = unit(heightPanel, "in"))
        
        ggsave(p_fixed, filename = file.path(folderToSavePlot,paste0('examplePlot_',peakSel, '_',cytoSel.short,'_covariateCorrectedSingleTp_flipped_',indexSize,'.pdf')), 
               width = widthSel, height=heightSel, device = cairo_pdf)
      }
    }
  }
}