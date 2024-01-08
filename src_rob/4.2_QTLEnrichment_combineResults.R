library(R.utils)
library(data.table)
library(stringr)
library(Cairo)
library(ggplot2)

###Get the variables we need
library(config)
Sys.setenv(R_CONFIG_ACTIVE = "300bcg")
config <- config::get()
for (x in names(config)){eval(parse(text=paste0(x,"='",config[[x]],"'")))}

source(file.path('300bcg_generalFxns.R'))

plotTopPathwaysNrSignif <- function(selPvals.df, topToPlot=20, titleSel=NULL,colorOutline='#937860'){
  #selPvals.df - df with one column that contains the pvalues of interest
  #colSel - selected column name with the pvalues of interest
  nrSignif = rowSums(selPvals.df<.05)
  minPvals = apply(selPvals.df, 1, min)
  nrSignif.df = cbind.data.frame('nrSignif'=nrSignif,'minPvals'=minPvals[names(nrSignif)])
  newOrder = with(nrSignif.df, order(-nrSignif, minPvals))
  nrSignif.df = nrSignif.df[newOrder,]
  nrSignif.df$pathways = rownames(nrSignif.df)
  
  ##Sort by two different values: first by number of signif and after that by minimal p-value
  #Select the top
  if (nrow(nrSignif.df)<topToPlot){
    topToPlot = nrow(nrSignif.df)
  }
  nrSignif.df$percSignif = 100 * ( nrSignif.df$nrSignif/ncol(selPvals.df) )
  nrSignifToPlot.df = nrSignif.df[c(1:topToPlot),]
  nrSignifToPlot.df$pathway = factor(nrSignifToPlot.df$pathways, levels = rev(nrSignifToPlot.df$pathways))
  
  maxVal = max(abs(nrSignifToPlot.df$percSignif))
  nrSignifToPlot.df$direction='Up'
  p1 = plotTopPathwaysCounts(nrSignifToPlot.df, pathwaySelName=titleSel, colorOutline=colorOutline)
  
  return(list('plot'=p1,'nrSignif.df'=nrSignif.df))
}

##Load all data
allRelevantDataList = loadLukasMetaDataPbmc(allMetaDataPath)
list2env(allRelevantDataList, .GlobalEnv)

allPathwayFiles = list.files(pascalOutputFolder) ##make sure this is set to the correct folder where pascal saves the output

#Define the groups we want to look at
cytokineGroups = defineSubgroupsCytokines()

#Define the combinations we want to look at, for both the fold-changes and the absolute values
combinationsToCheck = t(as.data.frame(list( c("cytoData.good.fc","V3"),
                                            c("cytoData.good.fc","V2"),
                                            c("cytoData.good.log2","V1"),
                                            c("circMarkerData.log2","V1")
)))

print(combinationsToCheck)
for (rowIndex in c(1:nrow(combinationsToCheck))){
  
  dataTypeSel=combinationsToCheck[[rowIndex,1]]
  timepointSel=combinationsToCheck[[rowIndex,2]]
  
  dataSel = allRelevantDataList[[as.character(dataTypeSel)]]
  selVariables = colnames(dataSel)
  
  #'KeggAndGo'=c("GO_Biological_Process_2018",'KEGG_2019_Human'),
  pathwaySetList = list('KEGG_2019_Human'=c('KEGG_2019_Human'), 
                        'GO_Biological_Process_2018' = c("GO_Biological_Process_2018") )
  
  #Create a combined list with all the pathway results
  combinedPathwayList = list()
  methodSel = 'max'
  combinedPathwayList = list()
  for (pathwaySelName in names(pathwaySetList)){
    ###Start combining the data for this pathway
    pathwaySetSel = pathwaySetList[[pathwaySelName]]
    
    print(methodSel)
    print(pathwaySetSel)
    combinedData = NULL
    for (varSel in selVariables){
      
      selData = NULL
      for (pathwaySel in pathwaySetSel){
        selFile = allPathwayFiles[grep(pattern = paste0('^',varSel,'_.*','\\_\\_',timepointSel,'\\_\\_',pathwaySel,"\\.",methodSel,'\\.txt$'),allPathwayFiles)]
        readData = fread(file.path(pascalOutputFolder,selFile))
        selData.part = readData[,c('Name','empPvalue')]
        colnames(selData.part) = c('Pathway',paste('empPvalue',varSel,sep='_'))
        if (is.null(selData)){
          selData = selData.part
        }else{
          selData = rbind.data.frame(selData,selData.part)
        }
      }
      
      if (is.null(combinedData)){
        combinedData = selData
      }else{
        combinedData = base::merge(combinedData, selData, by = "Pathway", all=T)
      }
    }
    combinedPathwayList[[pathwaySelName]] = combinedData
    
    ###Now start making the plots
    combinedData = as.data.frame(combinedData)
    rownames(combinedData) = combinedData[,1]
    combinedData[,1] <- NULL
    for (subGroupSel in names(cytokineGroups[[dataTypeSel]])){
      #Take the subgroup we are interested in
      print(subGroupSel)
      colsToGrep = unlist(lapply(cytokineGroups[[dataTypeSel]][[subGroupSel]], function(x) grep(pattern = paste0('.*\\_',x,'$') , colnames(combinedData))))
      combinedData.sel = combinedData[,colsToGrep,drop=F]
      nrSignif = rowSums(combinedData.sel<.05)
      combinedData.sel = combinedData.sel[order(nrSignif, decreasing = T),]
      rownames(combinedData.sel) = processPathwayNames(rownames(combinedData.sel))
      
      dirToSaveEnrichImages = file.path(image_dir_300BCG,'manuscriptImagesRaw','genetics','pathwayResults','pascal',paste(dataTypeSel,timepointSel,sep='_'))
      mkdirs(dirToSaveEnrichImages)
      for (topToPlot in c(5,10)){
        #Reformat data and make the plot
        plotList = plotTopPathwaysNrSignif(combinedData.sel, topToPlot=topToPlot, titleSel=pathwaySelName, colorOutline ='#937860')
        plot1 = plotList[['plot']]
        nrSignif.df = plotList[['nrSignif.df']]
        
        #Define panel and figure sizes
        if (dataTypeSel=="circMarkerData.log2"){
          widthSel=3.4
          widthPanel = .6
        }else{
          widthSel=4
          widthPanel = 1.2
        }
        
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
        #Set the panel size so that they are equal between both the genetics and epigenetics
        p_fixed <- set_panel_size(plot1,
                                  width  = unit(widthPanel, "in"),
                                  height = unit(heightPanel, "in"))
        
        print(dirToSaveEnrichImages)
        openxlsx::write.xlsx(nrSignif.df, file = file.path(dirToSaveEnrichImages,paste0(dataTypeSel,'_',subGroupSel,'_',pathwaySelName,'_',methodSel,'.xlsx')), rowNames=T)
        ggsave(p_fixed, filename = file.path(dirToSaveEnrichImages,paste0(dataTypeSel,'_',subGroupSel,'_',pathwaySelName,'_',methodSel,'_topNr=',topToPlot,'.pdf')),
               width=widthSel, height=heightSel)#, device = cairo_pdf)
      }
    }
  }
  
  ####Save this combined file as a supplement
  ##Combine the two pathway sets
  pathwayListBothLibs = NULL
  for (pathwaySelName in names(pathwaySetList)){
    selData = combinedPathwayList[[pathwaySelName]]
    selData.m = melt.data.table(selData, id.vars = "Pathway")
    selData.m$variable = gsub('^empPvalue_','',selData.m$variable)
    selData.m$variable = gsub('.fc$','',selData.m$variable)
    selData.m$variable = mapCytoNamesForSupplement(selData.m$variable,dataTypeSel)
    selData.m[[definitionsForSupplement[['library']]]] = pathwaySelName
    colnames(selData.m) = c(definitionsForSupplement[['gene_set']],definitionsForSupplement[[paste(dataTypeSel,'LONG',sep='_')]],definitionsForSupplement[['pval']],definitionsForSupplement[['library']])
    orderOfCols = c(definitionsForSupplement[['library']],definitionsForSupplement[['gene_set']],definitionsForSupplement[[paste(dataTypeSel,'LONG',sep='_')]],definitionsForSupplement[['pval']])
    selData.m = selData.m[,..orderOfCols]
    if (is.null(pathwayListBothLibs)){
      pathwayListBothLibs = selData.m
    }else{
      pathwayListBothLibs = rbind(pathwayListBothLibs,selData.m)
    }
  }
  
  pathwayListBothLibs = pathwayListBothLibs[order(pathwayListBothLibs[[definitionsForSupplement[['pval']]]]),]
  
  ## save xlsx
  mkdirs(file.path(suppTableDir,'xlsx_small','raw'))
  fileNameXlsx = file.path(suppTableDir,'xlsx_small','raw',paste0('TableS6_',tableAdditionToFileNames[['TableS6']],'.xlsx'))
  if (!(file.exists(fileNameXlsx))){
    wb <- openxlsx::createWorkbook()
  }else{
    wb <- openxlsx::loadWorkbook(fileNameXlsx)
  }
  sheetName = paste0(definitionsForSupplement[["SNP"]],'_',definitionsForSupplement[[dataTypeSel]],'_',definitionsForSupplement[[timepointSel]])
  print(sheetName)
  if (!(sheetName %in% wb$sheet_names)){
    openxlsx::addWorksheet(wb = wb, sheet = sheetName)
  }
  
  bold_style <- openxlsx::createStyle(textDecoration = "Bold")
  # make non-numeric column bold
  selCols = which(unlist(lapply(pathwayListBothLibs,class))!="numeric")
  openxlsx::addStyle(wb, style = bold_style, sheet = sheetName, rows = 1:(nrow(pathwayListBothLibs)+1), cols = selCols, gridExpand = TRUE)
  
  # create style, in this case bold header
  header_st <- openxlsx::createStyle(textDecoration = "Bold")
  
  openxlsx::writeData(wb = wb, sheet = sheetName, x = pathwayListBothLibs, colNames = T, rowNames = F, headerStyle = header_st)
  openxlsx::saveWorkbook(wb = wb, file = fileNameXlsx, overwrite = T)  
}