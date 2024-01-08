library("R.utils")
library('openxlsx')
library(ggplot2)
library(stringr)
library(data.table)
library(dplyr)
library(stats)
library(R.utils)
library(Cairo)
library(showtext)

print('starting')

###Get the variables we need
library(config)
Sys.setenv(R_CONFIG_ACTIVE = "300bcg")
config <- config::get()
for (x in names(config)){eval(parse(text=paste0(x,"='",config[[x]],"'")))}

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

source('300bcg_generalFxns.R')

###Example submissions
#sbatch 300bcg_atacGeneticsBayesian_plots.batch cytoData.good.log2 V1
#sbatch 300bcg_atacGeneticsBayesian_plots.batch circMarkerData.log2 V1
#sbatch 300bcg_atacGeneticsBayesian_plots.batch cytoData.good.fcCorr V3
#sbatch 300bcg_atacGeneticsBayesian_plots.batch cytoData.good.fcCorr V2
###Example args
#args = c('cytoData.good.fcCorr','V2')
#args = c('cytoData.good.fcCorr','V3')

args = base::commandArgs(trailingOnly=TRUE)
print(args)
stopifnot(length(args)==2)
dataTypeSel = args[[1]]
timepointSel = args[[2]]

##Load all metadata and data (the outcome factor)
allRelevantDataList = loadLukasMetaDataPbmc(allMetaDataPath)
list2env(allRelevantDataList, .GlobalEnv)

plotAllPercExplainedVar = T
if (plotAllPercExplainedVar){
  library(data.table) # to read the csv file
  
  processSigma <- function(sigmaValues){
    sigmaValues[sigmaValues>10]=0
    return(sigmaValues)
  }
  
  ####Settings
  tailNr = 5000
  subFolderSel = 'run1'
  folderToSaveImages = file.path(image_dir_300BCG,'manuscriptImagesRaw','atacSeqGeneticsCombined', 'percOfVariance_plots')
  mkdirs(folderToSaveImages)
  
  folderToSaveXlsx = file.path(suppTableDir,'xlsx_small','raw')
  mkdirs(folderToSaveXlsx)
  fileNameXlsx = file.path(folderToSaveXlsx,paste0('TableS7_',tableAdditionToFileNames[['TableS7']],'_pt1.xlsx'))
  
  nrOfVars = 50
  regressOut = T
  
  dataFrameSummary.list = list()
  
  for (geneticsInclusion in c('withGenetics','withoutGenetics')){
    dataFrameSummary = as.data.frame(matrix(NA, ncol=3,nrow=0))
    colnames(dataFrameSummary) = c('factor','sigma','FVE')
    print('genetics')
    
    for (seedSel in c(1,2,3,4,5)){
      #####SETTINGS
      topNrHits = as.numeric(nrOfVars)
      regressOutOfAtac = as.logical(regressOut)
      regressOutOfCells = as.logical(regressOut)
      regressOutOfY = F
      filterTopHits = T
      maxIter = 1000000
      burnIn = 900000
      geneticsPCutoff = "0"
      atacPCutoff = "0"
      
      regressOutStr = substr(as.character(regressOut), start=1, stop=1)
      argsSel = c(dataTypeSel,timepointSel,as.character(seedSel),as.character(as.integer(maxIter)),as.character(as.integer(burnIn)),
                  regressOutStr,regressOutStr,'F','T',as.character(nrOfVars),atacPCutoff,geneticsPCutoff)
      
      folderToSaveBayesFiles = file.path(data_dir_300BCG,'genetics','atacGeneticsCompared','bayesRR',subFolderSel,paste(dataTypeSel, timepointSel,sep='_'))
      
      #####LOAD
      for (selVarIndex in c(1:ncol(allRelevantDataList[[dataTypeSel]]))){
        selVar = colnames(allRelevantDataList[[dataTypeSel]])[[selVarIndex]]
        print(selVar)
        
        fileSel = list.files(path = folderToSaveBayesFiles, pattern = paste0(selVar,"_",paste(argsSel, collapse='_'),'_',selVarIndex,'_',geneticsInclusion,'\\.*'))
        if (length(fileSel)==0){next}
        bayesROutputFile = file.path(folderToSaveBayesFiles,fileSel)
        C<-fread(bayesROutputFile)
        
        selSigmas = names(C)[grep('^sigma.*', x = names(C))]
        for (sigmaSel in selSigmas){
          #Get values for selected sigma
          valuesSigmaSel = processSigma(C[[sigmaSel]])
          
          #Calculate sum of all sigmas
          valsAllSigmas = lapply(selSigmas, function(x) processSigma(C[[x]]))
          valsAllSigmas.df = as.data.frame(do.call(rbind, valsAllSigmas))
          summedSigmas = colSums(valsAllSigmas.df)
          
          #Calculate median and quartiles
          selVals.df = as.data.frame(tail(valuesSigmaSel,1000)/tail(summedSigmas,1000))
          colnames(selVals.df) = 'FVE'
          selVals.df$factor = selVar
          selVals.df$sigma = sigmaSel
          
          dataFrameSummary = rbind.data.frame(dataFrameSummary,selVals.df[,colnames(dataFrameSummary)])
        }
      }
    }
    dataFrameSummaryRaw = dataFrameSummary
    
    ##sigma mapping
    if (timepointSel=='V1'){
      sigmaMappingVector = c('Residuals','Host','Cell freq.','Chromatin','Genetics')
      names(sigmaMappingVector) = c("sigmaE","sigmaG[1]","sigmaG[2]","sigmaG[3]","sigmaG[4]")
    }else{
      sigmaMappingVector = c('Residuals','Host','Season','Cell freq.','Chromatin','Genetics')
      names(sigmaMappingVector) = c("sigmaE","sigmaG[1]","sigmaG[2]","sigmaG[3]","sigmaG[4]","sigmaG[5]")
    }
    
    mappedSigmas = unlist(lapply(dataFrameSummary$sigma, function(x) sigmaMappingVector[[x]]))
    dataFrameSummary$sigma = mappedSigmas
    dataFrameSummary$sigma = factor(dataFrameSummary$sigma, levels = c('Host','Season','Cell freq.','Genetics','Chromatin','Residuals'))
    oldLabs = as.character(dataFrameSummary$factor)
    oldLabs = gsub('\\.fc.*','',oldLabs)
    
    ##cytokine mapping
    if (startsWith(dataTypeSel,'cyto')){
      color.by = 'stimulus'
      newLabs = mapCytokineNamesForPlot(oldLabs, color.by=color.by)
      dataFrameSummary$factorMapped = newLabs 
      
      dataFrameSummary$factorMapped = factor(dataFrameSummary$factorMapped, levels = rev(orderToPlotCytokines))
      dataFrameSummary$stimulus = gsub(pattern = '.*\\(([A-Z,a-z]{2,3})\\)$',replacement = '\\1',dataFrameSummary$factorMapped)
      
      colorsSelList = c('#55a868',"#dd8452","#8172b3","#da8bc3")
      names(colorsSelList) = c("Mt","Sa","Ca","LPS")
      custom.col.signif=c()
      orderedStims = unique(dataFrameSummary$stimulus)
      for (stimSel in orderedStims){
        stimShort = names(colorsSelList)[[grep(pattern = paste0('^',substr(stimSel,1,1),'.*'), names(colorsSelList))]]
        custom.col.signif = c(custom.col.signif,colorsSelList[[stimShort]])
      }
    }else{
      color.by = 'other'
      newLabs = mapCytokineNamesForPlot(str_split_fixed(oldLabs,'_',2)[,1], color.by=color.by)
      dataFrameSummary$factorMapped = newLabs 
      dataFrameSummary$factorMapped = factor(dataFrameSummary$factorMapped, levels = rev(orderToPlotInflaMark))
      dataFrameSummary$stimulus = dataFrameSummary$factorMapped
      custom.col.signif = rep(c( "#4c72b0", "#769ad1"), length(unique(dataFrameSummary$sigma)))
    }
    
    #calculate median per factor
    library(dplyr)
    summarizedVals = dataFrameSummary%>%
      dplyr::group_by(factorMapped,sigma)%>% 
      dplyr::summarise(Mean=mean(FVE), Max=max(FVE), Min=min(FVE), Median=median(FVE), Std=sd(FVE))
    
    medianOverview = as.data.frame(summarizedVals[which(summarizedVals$sigma=='Chromatin'),c('factorMapped','Median')])
    
    #https://stackoverflow.com/questions/12570816/ggplot-scatter-plot-of-two-groups-with-superimposed-means-with-x-and-y-error-bar
    dataFrameSummary.meanSd = dataFrameSummary %>% group_by(factorMapped, sigma) %>% summarise_at(vars(FVE), list(min=min, Q1=~quantile(., probs = 0.25),
                                                                                                                  median=median, Q3=~quantile(., probs = 0.75),
                                                                                                                  max=max))
    dataFrameSummary.meanSd$stimulus = gsub(pattern = '.*\\(([A-Z,a-z]{2,3})\\)$',replacement = '\\1',dataFrameSummary.meanSd$factorMapped)
    
    fillList = list()
    fillList[['stimulus']] = "stimulus"
    fillList[['other']] = "sigma"
    
    textSize = 12
    textSizeLabs = 9
    
    dataFrameSummary.list[[geneticsInclusion]] = dataFrameSummary.meanSd

    ###Also make a general plots showing distribution of medians
    dataFrameViolin = as.data.frame(dataFrameSummary.meanSd)
    dp <- ggplot(dataFrameViolin, aes(x=sigma, y=median, fill=sigma)) + 
      geom_violin(trim=TRUE, scale='width')+
      ylab("Variance explained (%)")
    dp = dp +       theme_bw(base_family = "Helvetica") + 
      theme( 
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        legend.position = "none",
        text = element_text(family = "Helvetica"),
        axis.text = element_text(family = "Helvetica", size=textSize, colour = 'black'),
        axis.text.x = element_text(angle = 90,hjust=0.95,vjust=0.2),
        axis.title.y = element_text(size=textSize, color='black',family = "Helvetica"),
        strip.text = element_text(size=textSize, color='black',family = "Helvetica"),#added
        axis.title.x = element_blank()
      ) + scale_y_continuous(labels = scales::percent_format(accuracy = 1)) 
    
    
    heightSel = 3.5
    heightPanel =  1 
    widthSel=5
    widthPanel = 2.5
    library(egg)
    #Fix the dimensions of the panel, so all of them nicely align
    dp_fixed <- set_panel_size(dp,
                              width  = unit(widthPanel, "in"),
                              height = unit(heightPanel, "in"))
    
    ggsave(dp_fixed, filename = file.path(folderToSaveImages,paste0(paste(argsSel,collapse = '_'),'_',geneticsInclusion,'_boxplots_cairo.pdf')),
           width=widthSel,height=heightSel, units = "in")
    
    
    p = ggplot(dataFrameSummary.meanSd, aes_string(y = "factorMapped", x = "median", fill=fillList[[color.by]])) + 
      geom_bar(stat='identity') + 
      geom_errorbar(aes(xmin=Q1, xmax=Q3), width=.2,
                    position=position_dodge(.9))  +
      facet_wrap(~sigma, nrow=1)  +
      xlab("Variance explained (%)") + 
      ylab("") + 
      theme_bw(base_family = "Helvetica") + 
      theme( 
        panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank(),
        legend.position = "none",
        text = element_text(family = "Helvetica"),
        axis.text = element_text(family = "Helvetica", size=textSizeLabs, colour = 'black'),
        axis.title.x = element_text(size=textSize, color='black',family = "Helvetica"),
        strip.text = element_text(size=textSize, color='black',family = "Helvetica"),#added
        axis.title.y = element_blank()
      ) +
      scale_fill_manual(values = custom.col.signif) +
      scale_x_continuous(limits = c(0,.8), breaks = c(0,.5),labels = scales::percent_format(accuracy = 1)) #+
    
    if (dataTypeSel %in% c('circMarkerData.log2')){
      widthSel = 5.2412 * 2 * 0.393701
      heightSel = 14.0278 * 2 * 0.393701
    }else if (dataTypeSel %in% c("cytoData.good.fcCorr")){  
      widthSel = 5.2412 * 2 * 0.393701 * (7/6)
      heightSel = 6.7089 * 2 * 0.393701
    }else{
      widthSel = 5.2412 * 2 * 0.393701
      heightSel = 6.7089 * 2 * 0.393701
    }
  
    CairoFonts(
      regular="Helvetica:style=Regular",
      bold="Helvetica:style=Bold",
      italic="Helvetica:style=Italic"
    )

    ggsave(p, filename = file.path(folderToSaveImages,paste0(paste(argsSel,collapse = '_'),'_',geneticsInclusion,'_cairo.pdf')),
           width=widthSel,height=heightSel, units = "in", device = grDevices::cairo_pdf)
    
  }
  
  #####Save supplemental tables
  #First save the main supplement
  #Order: cytokine, factor, Q1, median, Q3
  selSaveSet = dataFrameSummary.list[['withGenetics']]
  selSaveSet = selSaveSet[,c("factorMapped","sigma","Q1","median","Q3")]
  colnames(selSaveSet) = c(definitionsForSupplement[[paste(dataTypeSel,'LONG',sep='_')]],definitionsForSupplement[["factor_group"]],
                           definitionsForSupplement[["FVE_Q1"]],definitionsForSupplement[["FVE_median"]],definitionsForSupplement[["FVE_Q3"]])
  selSaveSet[[definitionsForSupplement[["factor_group"]]]] = factor(selSaveSet[[definitionsForSupplement[["factor_group"]]]], levels = unique(selSaveSet[[definitionsForSupplement[["factor_group"]]]]))
  selSaveSet[[definitionsForSupplement[[paste(dataTypeSel,'LONG',sep='_')]]]] = factor(selSaveSet[[definitionsForSupplement[[paste(dataTypeSel,'LONG',sep='_')]]]], 
                                                                                       levels = rev(levels(selSaveSet[[definitionsForSupplement[[paste(dataTypeSel,'LONG',sep='_')]]]])))
  
  selSaveSet = selSaveSet[order(unlist(selSaveSet[,definitionsForSupplement[[paste(dataTypeSel,'LONG',sep='_')]]]),unlist(selSaveSet[,definitionsForSupplement[["factor_group"]]])),]
  
  sheetName = paste('varPart',dataTypeNameMapSupplement[[dataTypeSel]],mapVisitShortNameToFinal[[timepointSel]],sep='_')
  if (endsWith(suffix = '.fcCorr', dataTypeSel)){
    sheetName = paste0(sheetName,'FC')
  }
  
  if (!(file.exists(fileNameXlsx))){
    wb <- openxlsx::createWorkbook()
  }else{
    wb <- openxlsx::loadWorkbook(fileNameXlsx)
  }
  
  if (!(sheetName %in% wb$sheet_names)){
    openxlsx::addWorksheet(wb = wb, sheet = sheetName)
  }
  
  bold_style <- createStyle(textDecoration = "Bold")
  # make non-numeric column bold
  selCols = which(unlist(lapply(selSaveSet,class))!="numeric")
  addStyle(wb, style = bold_style, rows = 1:(nrow(selSaveSet)+1), sheet = sheetName, cols = selCols, gridExpand = TRUE)
  
  # create style, in this case bold header
  header_st <- openxlsx::createStyle(textDecoration = "Bold")
  
  openxlsx::writeData(wb = wb, sheet = sheetName, x = selSaveSet, colNames = T, rowNames = F, headerStyle = header_st)
  openxlsx::saveWorkbook(wb = wb, file = fileNameXlsx, overwrite = T)
  
  #####Merge and plot with and without genetics
  res.part1 = dataFrameSummary.list[['withGenetics']]
  res.part1$genetics_used = 'withGenetics'
  res.part2 = dataFrameSummary.list[['withoutGenetics']]
  res.part2$genetics_used = 'withoutGenetics'
  
  res.combined = rbind.data.frame(res.part1,res.part2)
  res.combined.sel= res.combined[res.combined$sigma=='Chromatin',]
  
  p = ggplot(res.combined.sel, aes_string(y = "factorMapped", x = "median", fill='genetics_used')) + 
    geom_bar(stat='identity',position="dodge") + 
    geom_errorbar(aes(xmin=Q1, xmax=Q3), width=.2,
                  position=position_dodge(.9))  +
    facet_wrap(~sigma, nrow=1)  +
    xlab("Variance explained (%)") + 
    ylab("") + 
    theme_bw(base_family = "Helvetica") + 
    theme( 
      panel.grid.major.y = element_blank(),
      panel.grid.minor.y = element_blank(),
      axis.text = element_text(family = "Helvetica", size=textSizeLabs, colour = 'black'),
      axis.title.x = element_text(size=textSize, color='black',family = "Helvetica"),
      axis.title.y = element_blank()
    ) +
    scale_fill_manual(name = "Subset", labels = c("Including genetics", "Excluding genetics"), values = c( "#4c72b0", "#769ad1")) +
    scale_x_continuous(limits = c(0,.8), breaks = c(0,.5),labels = scales::percent_format(accuracy = 1)) 
  
  p
  
  widthSel=5
  widthPanel = 2*5.1216*0.393701 *.50
  if (dataTypeSel=="circMarkerData.log2"){
    heightSel = 14.0278 * 2 * 0.393701
    heightPanel = heightSel*.92
  }else{
    #Set the panel size for this
    heightSel = 6.7089 * 2 * 0.393701
    heightPanel = 4.47
  }
  
  library(egg)
  #Fix the dimensions of the panel, so all of them nicely align
  p_fixed <- egg::set_panel_size(p,
                                 width  = unit(widthPanel, "in"),
                                 height = unit(heightPanel, "in"))
  
  CairoFonts(regular="Helvetica:style=Regular",
             bold="Helvetica:style=Bold",
             italic="Helvetica:style=Italic",
             bolditalic="Helvetica:style=Bold Italic,BoldItalic")
  
  ggsave(p_fixed, filename = file.path(folderToSaveImages,paste0(paste(argsSel,collapse = '_'),'_withAndWithoutGenetics.pdf')),
         width=widthSel,height=heightSel, units = "in", device = cairo_pdf)
  
  ##Save table
  selSaveSet = res.combined.sel
  selSaveSet = selSaveSet[,c("factorMapped","sigma","genetics_used","Q1","median","Q3")]
  colnames(selSaveSet) = c(definitionsForSupplement[[paste(dataTypeSel,'LONG',sep='_')]],
                           definitionsForSupplement[["factor_group"]],
                           "genetics_used",
                           definitionsForSupplement[["FVE_Q1"]],definitionsForSupplement[["FVE_median"]],definitionsForSupplement[["FVE_Q3"]])
  
  ##In the end save just the one without genetics
  selSaveSet = selSaveSet[which(selSaveSet[['genetics_used']]=='withoutGenetics'),]
  selSaveSet = selSaveSet[,c(definitionsForSupplement[[paste(dataTypeSel,'LONG',sep='_')]],
                             definitionsForSupplement[["factor_group"]],
                             definitionsForSupplement[["FVE_Q1"]],definitionsForSupplement[["FVE_median"]],definitionsForSupplement[["FVE_Q3"]])]
  #selSaveSet[["genetics_used"]] = factor(selSaveSet[["genetics_used"]], levels = unique(selSaveSet[["genetics_used"]]))
  selSaveSet[[definitionsForSupplement[[paste(dataTypeSel,'LONG',sep='_')]]]] = factor(selSaveSet[[definitionsForSupplement[[paste(dataTypeSel,'LONG',sep='_')]]]], 
                                                                                       levels = rev(levels(selSaveSet[[definitionsForSupplement[[paste(dataTypeSel,'LONG',sep='_')]]]])))
  selSaveSet = selSaveSet[order(unlist(selSaveSet[,definitionsForSupplement[[paste(dataTypeSel,'LONG',sep='_')]]])),]
  
  sheetNameSel = gsub('varPart','varPartExclSNP',sheetName)
  if (!(file.exists(fileNameXlsx))){
    wb <- openxlsx::createWorkbook()
  }else{
    wb <- openxlsx::loadWorkbook(fileNameXlsx)
  }
  
  if (!(sheetNameSel %in% wb$sheet_names)){
    openxlsx::addWorksheet(wb = wb, sheet = sheetNameSel)
  }
  
  bold_style <- createStyle(textDecoration = "Bold")
  # make non-numeric column bold
  selCols = which(unlist(lapply(selSaveSet,class))!="numeric")
  addStyle(wb, style = bold_style, rows = 1:(nrow(selSaveSet)+1), sheet = sheetNameSel, cols = selCols, gridExpand = TRUE)
  
  # create style, in this case bold header
  header_st <- openxlsx::createStyle(textDecoration = "Bold")
  
  openxlsx::writeData(wb = wb, sheet = sheetNameSel, x = selSaveSet, colNames = T, rowNames = F, headerStyle = header_st)
  openxlsx::saveWorkbook(wb = wb, file = fileNameXlsx, overwrite = T)
}