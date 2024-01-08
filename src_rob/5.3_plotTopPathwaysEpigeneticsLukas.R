library(data.table)
library(dplyr)
library(tidyr)
library("R.utils")
library('openxlsx')
library(ggplot2)
library(stringr)
library(data.table)
library(dplyr)
library(RColorBrewer)

###Get the variables we need
library(config)
Sys.setenv(R_CONFIG_ACTIVE = "300bcg")
config <- config::get()
for (x in names(config)){eval(parse(text=paste0(x,"='",config[[x]],"'")))}

source('plottingFxnsGeneral.R')
source('300bcg_generalFxns.R')

#Load all data
allRelevantDataList = loadLukasMetaDataPbmc(allMetaDataPath)
list2env(allRelevantDataList, .GlobalEnv)

#Define where to save the combined pathway data
folderToSave = file.path(data_dir_300BCG,'atacSeq','pathwayResults')
mkdirs(folderToSave)

#Define the different subsets of cytokines to look at
cytokineGroups = defineSubgroupsCytokines()

###Here I combine the pathway data per cytokine into complete sets of interest
reExtractPathways = T
if (reExtractPathways){
  ##We want to look at KEGG and GO
  for (pathwaySel in c("KEGG_2019_Human_min15_max500","GO_Biological_Process_2018_min15_max500")){
    ##Look at two different distance metrics that were used to define if a peak belongs to a gene
    distanceSel='GENE_AND_DISTAL_10kb'
    
    combinationsToCheck = t(as.data.frame(list(
      c("circMarkerData.log2","V1"),
      c("cytoData.good.log2","V1"),
      c("cytoData.good.fc","V2"),
      c("cytoData.good.fc","V3"))))
    
    print(combinationsToCheck)
    
    #Iterate over the different combinations of dataset and timepoint
    for (rowIndex in c(1:nrow(combinationsToCheck))){
      
      dataTypeSel=combinationsToCheck[[rowIndex,1]]
      timepointSel=combinationsToCheck[[rowIndex,2]]
      
      #Select the dataset
      selDataSet = allRelevantDataList[[dataTypeSel]]
      selVariables = colnames(selDataSet)
      
      #Define if we are looking at cytokines or circulating mediators
      if (startsWith(x = dataTypeSel, prefix = 'cytoData')){
        selData = 'cytokines'
      }else{
        selData = 'CM'
      }
      
      #Check if we are looking at fold-changes or absolute values
      if (endsWith(x = dataTypeSel, suffix = '.fc')){
        selFolder = file.path(folderWithTrainingAtacResults,selData)
      }else{
        selFolder = file.path(folderWithBaselineAtacResults,selData)
      }
      
      ##Go over each subfolder with the Enrichr data
      fullEnrichrSetUp = NULL
      fullEnrichrSetDown = NULL
      for (varSel in colnames(selDataSet)){
        print(varSel)
        #Check the relevant folder
        folderWithInfo = list.files(path = file.path(selFolder), pattern = paste0('(\\_|\\.)',timepointSel,'(\\_|\\.).*',gsub('\\.fc$','',varSel),'(\\_good)?$'), include.dirs = T)
        stopifnot(length(folderWithInfo)==1)
        
        #Enrichr folder
        folderWithInfoPerVar = file.path(selFolder, folderWithInfo, "Enrichr")
        
        #Load the data -- we are looking at the pathway results for enrichr, where we defined "signficant" in the fisher test as the top 1000 genes
        fileSelUp = list.files(path = folderWithInfoPerVar, pattern = paste0("^enrichr\\_results.*p\\.value\\.top\\_1000\\.up\\.",distanceSel,"\\.tsv$"))
        fileSelDown = list.files(path = folderWithInfoPerVar, pattern = paste0("^enrichr\\_results.*p\\.value\\.top\\_1000\\.down\\.",distanceSel,"\\.tsv$"))
        stopifnot(length(fileSelUp)==length(fileSelDown))
        stopifnot(length(fileSelUp)==1)
        #Load the data - separetely for up and down-regulated genes
        dataUp = fread(file.path(folderWithInfoPerVar,fileSelUp))
        dataUp = dataUp[which(dataUp$Gene_set==pathwaySel),c("Term","P-value")]
        colnames(dataUp) = c("Term",varSel)
        
        dataDown = fread(file.path(folderWithInfoPerVar,fileSelDown))
        dataDown = dataDown[which(dataDown$Gene_set==pathwaySel),c("Term","P-value")]
        colnames(dataDown) = c("Term",varSel)
        
        #combine
        if (is.null(fullEnrichrSetUp)){
          fullEnrichrSetUp = dataUp
          fullEnrichrSetDown = dataDown
        }else{
          fullEnrichrSetUp = base::merge(fullEnrichrSetUp,dataUp,by="Term",all=T)
          fullEnrichrSetDown = base::merge(fullEnrichrSetDown,dataDown,by="Term",all=T)
        }
      }
      
      dataTypeSel=combinationsToCheck[[rowIndex,1]]
      timepointSel=combinationsToCheck[[rowIndex,2]]
      
      #Save the full set
      fwrite(fullEnrichrSetDown, file = file.path(folderToSave,paste0(dataTypeSel,'_',mapVisitShortNameToFinal[[timepointSel]],'_',pathwaySetNameMapShort[[pathwaySel]],'_',distanceNameMapShort[[distanceSel]],'_pathwaysLukas_Down.csv')))
      fwrite(fullEnrichrSetUp, file = file.path(folderToSave,paste0(dataTypeSel,'_',mapVisitShortNameToFinal[[timepointSel]],'_',pathwaySetNameMapShort[[pathwaySel]],'_',distanceNameMapShort[[distanceSel]],'_pathwaysLukas_Up.csv')))
      
      saveFullSupplement = F
      if (saveFullSupplement){
        #############Supplement################
        ##Change the names of the circulating markers to match the names we use in the plots
        
        cytokineSel.up = gsub('\\.fc$','',colnames(fullEnrichrSetUp))[2:ncol(fullEnrichrSetUp)]
        colnames(fullEnrichrSetUp)[2:ncol(fullEnrichrSetUp)] = mapCytoNamesForSupplement(cytokineSel.up,dataTypeSel)
        
        cytokineSel.down = gsub('\\.fc$','',colnames(fullEnrichrSetDown))[2:ncol(fullEnrichrSetDown)]
        colnames(fullEnrichrSetDown)[2:ncol(fullEnrichrSetDown)] = mapCytoNamesForSupplement(cytokineSel.down,dataTypeSel)
        
        #Save to the same set to the path I use for saving the other supplementary tables
        #replace the "missing" values with p-value 1. These were not calculated by Lukas, since these had no significant hits, therefore no fisher results
        folderToSaveCsv = file.path(image_dir_300BCG,'manuscriptSupplementTables','atacSeq','pathwayResultsCsv',paste(dataTypeSel,timepointSel,sep='_'))
        mkdirs(folderToSaveCsv)
        fullEnrichrSetDown.filled = fullEnrichrSetDown
        fullEnrichrSetDown.filled[is.na(fullEnrichrSetDown)] = 1
        fullEnrichrSetUp.filled = fullEnrichrSetUp
        fullEnrichrSetUp.filled[is.na(fullEnrichrSetUp.filled)] = 1
        fwrite(fullEnrichrSetDown.filled, file = file.path(folderToSaveCsv,paste0(dataTypeSel,'_',mapVisitShortNameToFinal[[timepointSel]],'_',pathwaySetNameMapShort[[pathwaySel]],'_',distanceNameMapShort[[distanceSel]],'_pathwaysAtac_Down.csv')))
        fwrite(fullEnrichrSetUp.filled, file = file.path(folderToSaveCsv,paste0(dataTypeSel,'_',mapVisitShortNameToFinal[[timepointSel]],'_',pathwaySetNameMapShort[[pathwaySel]],'_',distanceNameMapShort[[distanceSel]],'_pathwaysAtac_Up.csv')))
        
        ##Also save xlsx
        folderToSaveXlsx = file.path(image_dir_300BCG,'manuscriptSupplementTables','atacSeq','pathwayResultsXlsx')
        mkdirs(folderToSaveXlsx)
        
        fileNameXlsx = file.path(folderToSaveXlsx,'atacPathwayResults.xlsx')
        if (!(file.exists(fileNameXlsx))){
          wb <- openxlsx::createWorkbook()
        }else{
          wb <- openxlsx::loadWorkbook(fileNameXlsx)
          
        }
        
        sheetNameUp = paste0(dataTypeNameMapShort[[dataTypeSel]],'_',mapVisitShortNameToFinal[[timepointSel]],'_',pathwaySetNameMapShort[[pathwaySel]],'_',distanceNameMapShort[[distanceSel]],'_up')
        sheetNameDown = paste0(dataTypeNameMapShort[[dataTypeSel]],'_',mapVisitShortNameToFinal[[timepointSel]],'_',pathwaySetNameMapShort[[pathwaySel]],'_',distanceNameMapShort[[distanceSel]],'_down')
        print(sheetNameUp)
        print(sheetNameDown)
        
        if (!(sheetNameUp %in% wb$sheet_names)){
          openxlsx::addWorksheet(wb = wb, sheet = sheetNameUp)
          openxlsx::addWorksheet(wb = wb, sheet = sheetNameDown)
        }
        
        openxlsx::writeData(wb = wb, sheet = sheetNameUp, x = fullEnrichrSetUp.filled, colNames = T, rowNames = F)
        openxlsx::writeData(wb = wb, sheet = sheetNameDown, x = fullEnrichrSetDown.filled, colNames = T, rowNames = F)
        openxlsx::saveWorkbook(wb = wb, file = fileNameXlsx, overwrite = T)
      }
    }
  }
}

###We are going to plot, just KEGG, just GO, and the combination of KEGG+GO
pathwaySetList = list('KeggAndGo'=c("Go",'Kegg'),
                      'Kegg'=c('Kegg'), 
                      'Go' = c("Go") )

##Make the plots
plotTopPathways = T
if (plotTopPathways){
  print('plotTopPathways')
  
  for (topToPlot in c(5,10)){
    for (pathwaySelName in names(pathwaySetList)){
      pathwaySetSel = pathwaySetList[[pathwaySelName]]
      
      #Again, check the two different distance metrics for mapping regions to genes
      distanceSel='GENE_AND_DISTAL_10kb'
      combinationsToCheck = t(as.data.frame(list(
        c("cytoData.good.log2","V1"),
        c("cytoData.good.fc","V3"),
        c("cytoData.good.fc","V2"),
        c("circMarkerData.log2","V1")
      )))
      print(combinationsToCheck)
      for (rowIndex in c(1:nrow(combinationsToCheck))){
        
        dataTypeSel=combinationsToCheck[[rowIndex,1]]
        timepointSel=combinationsToCheck[[rowIndex,2]]
        
        selDataSet = allRelevantDataList[[dataTypeSel]]
        selVariables = colnames(selDataSet)
        
        #Cytokines vs CM
        if (startsWith(x = dataTypeSel, prefix = 'cytoData')){
          selData = 'cytokines'
        }else{
          selData = 'CM'
        }
        
        if (endsWith(x = dataTypeSel, suffix = '.fc')){
          selFolder = file.path(folderWithTrainingAtacResults,selData)
        }else{
          selFolder = file.path(folderWithBaselineAtacResults,selData)
        }
        
        #Load all results into a list
        fullEnrichrSetList = list()
        for (directionSel in c('Down','Up')){
          fullEnrichrSetList[[directionSel]] = NULL
          for (pathwaySel in pathwaySetSel){
            
            fullEnrichrSet = fread( file.path(folderToSave,paste0(dataTypeSel,'_',mapVisitShortNameToFinal[[timepointSel]],'_',pathwaySel,'_',distanceNameMapShort[[distanceSel]],'_pathwaysLukas_',directionSel,'.csv')), data.table = F)
            rownames(fullEnrichrSet) = fullEnrichrSet$Term
            fullEnrichrSet[is.na(fullEnrichrSet)] = 1
            fullEnrichrSet$Term <- NULL
            if (is.null(fullEnrichrSetList[[directionSel]])){
              fullEnrichrSetList[[directionSel]] = fullEnrichrSet
            } else{
              fullEnrichrSetList[[directionSel]] = rbind.data.frame(fullEnrichrSetList[[directionSel]],fullEnrichrSet)
            }
          }
        }
        
        ##Now for the different subgroups we defined, merge the data and make the plots
        for (subGroupSel in names(cytokineGroups[[dataTypeSel]])){
          percentageSignif.list=list()
          allPvalDf = NULL
          for (directionSel in c('Down','Up')){
            fullEnrichrSet = fullEnrichrSetList[[directionSel]]
            print(subGroupSel)
            
            fullEnrichrSet.sel = fullEnrichrSet[,cytokineGroups[[dataTypeSel]][[subGroupSel]],drop=F]
            percSignif = rowSums(fullEnrichrSet.sel<=.05)/ncol(fullEnrichrSet.sel)
            names(percSignif) = rownames(fullEnrichrSet.sel)
            percSignif.df =as.data.frame(percSignif)
            colnames(percSignif.df) = directionSel
            percSignif.df$pathway = rownames(percSignif.df)
            percentageSignif.list[[directionSel]] = percSignif.df
            
            toMerge = fullEnrichrSet.sel
            toMerge$rn = rownames(toMerge)
            if (is.null(allPvalDf)){
              allPvalDf = toMerge
            }else{
              allPvalDf = base::merge(allPvalDf,toMerge,by='rn')
            }
          }
          rownames(allPvalDf) = allPvalDf$rn
          allPvalDf$rn <- NULL
          
          #Merge up and down
          mergedResults = base::merge(percentageSignif.list[['Up']], percentageSignif.list[['Down']], by='pathway',all=T)
          
          #Map to simpler names for the pathways
          rownames(mergedResults) =  mergedResults$pathway  
          mergedResults$pathway = processPathwayNames(mergedResults$pathway)
          
          #Check the number of significant (max between up and down)
          rowMaxVals = apply(mergedResults[,names(percentageSignif.list)], 1, max)
          #Sort first by the number of significant and second by the minimum p-value over all items
          minPvals = apply(allPvalDf, 1, min)
          minPvals = minPvals[names(rowMaxVals)]
          boundStats = cbind.data.frame('rowMaxVals'=rowMaxVals,'minPvals'=minPvals)
          newOrder = with(boundStats, order(-rowMaxVals, minPvals))
          
          mergedResults.order = mergedResults[newOrder,]
          mergedResults.order$Down = -mergedResults.order$Down
          mergedResults.order.sel = mergedResults.order[c(1:topToPlot),]
          head(mergedResults)
          
          #Reformat the df for the pot
          mergedResults.order.sel.m = melt(data = mergedResults.order.sel, id.vars = c("pathway"), measure.vars = c("Up", "Down"))
          mergedResults.order.sel.m$pathway = factor(mergedResults.order.sel.m$pathway, levels=rev(mergedResults.order[c(1:topToPlot),'pathway']))
          mergedResults.order.sel.m$percSignif = mergedResults.order.sel.m$value * 100
          mergedResults.order.sel.m$direction = mergedResults.order.sel.m$variable
          #Now make the plot
          p = plotTopPathwaysCounts(percentageResults.ordered=mergedResults.order.sel.m, pathwaySelName=pathwaySelName, colorOutline ='#ccb974')
          
          #Save the results
          dirToSaveEnrichImages = file.path(image_dir_300BCG,'manuscriptImagesRaw','atacSeq','pathwayResults',paste(dataTypeSel,timepointSel,sep='_'))
          mkdirs(dirToSaveEnrichImages)
          
          widthSel=4
          widthPanel = 1.2
          if (topToPlot==20){
            heightSel = 5
            heightPanel = 4
          }else if (topToPlot==10){
            heightSel = 2.3
            heightPanel = 2
          }else if (topToPlot==5){
            heightSel = 1.2
            heightPanel = 1
          }
          
          library(egg)
          #Fix the dimensions of the panel, so all of them nicely align
          p_fixed <- set_panel_size(p,
                                    width  = unit(widthPanel, "in"),
                                    height = unit(heightPanel, "in"))
          
          ggsave(p_fixed, filename = file.path(dirToSaveEnrichImages,paste0(dataTypeSel,'_',timepointSel,'_',subGroupSel,'_',pathwaySelName,'_',distanceSel,'_topNr=',topToPlot,'_lukas.pdf')),
                 width=widthSel, height=heightSel)
        }
      }
    }
  }
}
