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
source('plottingFxnsGeneral.R')
source('300bcg_generalFxns.R')

print('starting')

plotBoxesWithTopHits <- function(resultDf,   fdrCutOffs = c(.05,.1,.2,.3), minNrToAnnotate=15, 
                                 color.by='stimulus', annotateTopPerCyto=F, textSize=7,colorOutline ='#ccb974',rasterizePlot=T){
  #Create the "manhattan"-style boxplot showing the top results for the ATAC-seq analysis
  #*resultDf - Dataframe with the data to plot ordered by increasing p-value (from more to less significant)
  # Required columns:
  # -'P.Value'
  # -'adj.P.Val.total'
  # -'logPVal'
  # -'variable'
  # -'namesToPlot'
  #*fdrCutOffs - FDR cutoffs to show. All hits below the higest FDR included here ill be shown
  #*minNrToAnnotate - the mimimum number of hits to annotate with the gene closest to the region
  #*annotateTopPerCyto - if "T" in addition to the momimum number of annotated hits, per cytokine at least one hit will be annotated
  #*textSize - size of the annotations
  #*colorOutline - color of the outline around the plot
  
  library(ggrepel)
  #Check at which Cut-offs each region is significant
  fdrCutOffPlot = list()
  for (cutOffSel in fdrCutOffs){
    firstSignif = rev(which(resultDf$adj.P.Val.total<=cutOffSel))[[1]]
    nextNotSignif = firstSignif+1
    fdrCutOffPlot[[as.character(cutOffSel)]] = mean(unlist(resultDf[c(firstSignif,nextNotSignif),'p.value']))
  }
  
  #Now only plot the ones below the decided DR-value threshold
  resultDf = as.data.frame(resultDf)
  resultDf$jitter = NA
  resultDf[resultDf$p.value<=max(unlist(fdrCutOffPlot)),'jitter'] = 
    resultDf[resultDf$p.value<=max(unlist(fdrCutOffPlot)),'logPVal']
  
  #Mark the top X
  #Check how many adj.p are below .1
  #Annotate all hits below FDR .1 regardless of the minimum sest
  numberBelowPointOne = sum(resultDf$adj.P.Val.total<=.1)
  nrToAnnotate = max(minNrToAnnotate,numberBelowPointOne)
  
  #Do not annotate anything above the max fdrCutOffs
  numberBelowMaxFdr = sum(resultDf$adj.P.Val.total<=max(fdrCutOffs))
  nrToAnnotate = min(nrToAnnotate,numberBelowMaxFdr)
  #Unless there are litteraly 0
  if (numberBelowMaxFdr==0){
    nrToAnnotate =  minNrToAnnotate
  }
  
  #Select the top to annotate
  resultDf$topX = NA
  tpXSel = c(1:nrToAnnotate)
  resultDf[tpXSel,'topX'] = 
    resultDf[tpXSel,'logPVal']
  
  #If we want to also annotate the top hits per cytokine
  if (annotateTopPerCyto){
    nrToAnnotate=1
    for (cytoSel in unique(resultDf$variable)){
      rowsSel = which(resultDf$variable==cytoSel)
      if (nrToAnnotate<length(rowsSel)){
        rowsSel = rowsSel[c(1:nrToAnnotate)]
      }
      resultDf[rowsSel,'topX'] =resultDf[rowsSel,'logPVal']
    }
    addToName=paste0('top',nrToAnnotate,'AndTop',nrToAnnotate,'PerCyto')
  }else{
    addToName=paste0('top',nrToAnnotate)
  }
  
  #Keep just the variables with significant results
  resultDf = resultDf[resultDf$variable %in% unique(resultDf[!is.na(resultDf$jitter),'variable']),]
  
  #helvetica
  #Color palette
  if (color.by=='stimulus'){
    #For the cytokines we color by stimulus
    colorsSelList = c('#55a868',"#dd8452","#8172b3","#da8bc3")
    names(colorsSelList) = c("Mt","Sa","Ca","LPS")
    custom.col.signif=c()
    orderedStims = sort(unique(resultDf$variable))
    for (stimSel in orderedStims){
      stimShort = names(colorsSelList)[[grep(pattern = paste0('^',substr(stimSel,1,1),'.*'), names(colorsSelList))]]
      custom.col.signif = c(custom.col.signif,colorsSelList[[stimShort]])
    }
    custom.col.grey =  custom.col.signif
    custom.col.grey[custom.col.signif==colorsSelList[[1]]] = "#8c8c8c" 
    custom.col.grey[custom.col.signif==colorsSelList[[2]]] = "#b7b7b7" 
    custom.col.grey[custom.col.signif==colorsSelList[[3]]] = "#8c8c8c" 
    custom.col.grey[custom.col.signif==colorsSelList[[4]]] = "#b7b7b7" 
  }else{
    #For the circulating markers we color by bluescales
    custom.col.signif <- rep(c( "#4c72b0", "#769ad1"), length(unique(resultDf$variable)))
    custom.col.grey = rep(c( "#8c8c8c", "#b7b7b7"), length(unique(resultDf$variable)))
  }
  
  #Check if the association is positive or negative
  resultDf$associationDirection = NA
  resultDf$associationDirection[sign(resultDf$Coef)==1]='pos'
  resultDf$associationDirection[sign(resultDf$Coef)==-1]='neg'
  
  # Define which are outliers and inliers
  # From http://datacornering.com/how-to-create-boxplot-in-r-and-extract-outliers/
  resultDf$outlier <- ifelse(resultDf$logPVal %in% boxplot(resultDf$logPVal ~ resultDf$variable, plot = F)$out, resultDf$logPVal, NA)
  resultDf$inlier <- ifelse(is.na(resultDf$outlier), resultDf$logPVal, NA)
  #In the end we just want to plot those that are are not above any signficance line in grey
  resultDf$outlierToPlot = resultDf$outlier
  resultDf[which(resultDf$adj.P.Val.total<=max(fdrCutOffs)),'outlierToPlot'] = NA
  
  
  ####Start the plot
  #Jitter top hits slightly left and right
  pos <- position_jitter(width = 0.25,height = 0, seed=10)
  p <- ggplot(resultDf)
  p = p + geom_boxplot(outlier.shape = NA, aes(x = variable, y = logPVal, fill = variable))
  
  breaks = seq(0,to = ceiling(max(resultDf[,'logPVal'])/2 )*2, by=2)
  labels = breaks
  #Add lines for the FDR thresholds
  for (fdrCutOffName in names(fdrCutOffPlot)){
    logVal = -log10(fdrCutOffPlot[[fdrCutOffName]])
    breaks = c(breaks,logVal)
    labels = c(labels,fdrCutOffName)
    p = p + geom_hline(yintercept = logVal, linetype= "dashed", color = "#4c72b0")
  }
  
  #Plot the grey outling points
  set.seed(10)
  p = p + changePlotRasterization(geom_jitter(aes(x=variable, y=outlierToPlot), width = .25,
                      height = 0, shape=16, size=1.8, color = "#8c8c8c", alpha=.75))
  
  #Plot the points
  p = p + geom_point(position = pos, aes(x=variable, y=jitter, color=variable, fill=variable, shape=associationDirection), size=2.3, alpha=.75) + 
    scale_fill_manual(values = custom.col.signif) +
    scale_shape_manual(values = c("pos" = 24, "neg" = 25))

  #Plot the labels
  p = p + geom_label_repel(mapping = aes(x=variable, y=topX, color=variable, label = namesToPlotItalic), 
                           position = pos, min.segment.length=.1, size = textSize * (5/14.2256), parse=T
                           ) +  
    scale_color_manual(values = custom.col.signif)
  
  #Change the theme
  p = p + theme_bw(base_family = "Helvetica") + 
    theme(
      text = element_text(family = "Helvetica", size=textSize, colour = 'black'),
      axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, colour = 'black'), 
      legend.position = "none", 
      axis.text = element_text(size=textSize, color = 'black',family = "Helvetica"),
      axis.title.y = element_text(size=textSize, color='black',family = "Helvetica"),
      axis.title.x=element_blank()
    )
  
  p = p + ylab(expression(-log[10]*italic(P)*textstyle("-")*value))

  p = p + scale_y_continuous(breaks = breaks, labels = labels, expand = expansion(mult = c(0,0.05))) 
  
  #Some options that are different for the stimulated cytokines vs the circulating markers of inflammation
  oldLabs = ggplot_build(p)$layout$panel_params[[1]]$x$get_labels()
  print(paste('old labels = ', oldLabs))
  if (color.by=='stimulus'){
    cytokinesSel = str_split_fixed(oldLabs, pattern = '_', n = 5)[,4]
    mappedCytokines = unlist(lapply(cytokinesSel, function(x) mappingListCyto[[x]]))
    stimsShortMapped = mapLongCytoToShortCytoName(str_split_fixed(oldLabs, pattern = '_', n = 5)[,1])
    newLabs = paste0( mappedCytokines, ' (', stimsShortMapped, ')')
    #Reorder
    timePoint.fact = factor(str_split_fixed(oldLabs, pattern = '_', n = 5)[,2],levels = c('24h','7d'))
    stimsShortMapped.fact = factor(stimsShortMapped,levels=c('Ca','LPS','Sa','Mt'))
    mappedCytokines.fact = factor(mappedCytokines,sort(unique(mappedCytokines)))
    orderOfLabels = order(timePoint.fact, stimsShortMapped.fact,mappedCytokines.fact)
    
    print(paste0('limits = ',oldLabs[orderOfLabels]))
    print(paste0('newlabs = ',newLabs[orderOfLabels]))
    p = p + scale_x_discrete(limits = oldLabs[orderOfLabels], labels= newLabs[orderOfLabels])
  }else {
    mappedMarkers = unlist(lapply(oldLabs, function(x) mappingListCircMed[[x]]))
    orderedLabels = factor(mappedMarkers, levels = orderToPlotInflaMark)
    orderOfLabels = order(orderedLabels)
    p = p + scale_x_discrete(limits = oldLabs[orderOfLabels], labels= mappedMarkers[orderOfLabels])
  }
  
  return(list(p=p,addToName=addToName))
}


###Get the variables we need
library(config)
Sys.setenv(R_CONFIG_ACTIVE = "300bcg")
config <- config::get()
for (x in names(config)){eval(parse(text=paste0(x,"='",config[[x]],"'")))}

#Load cytokine and metadata
allRelevantDataList = loadLukasMetaDataPbmc(allMetaDataPath)
list2env(allRelevantDataList, .GlobalEnv)

##Load files to map stimuli and cytokines to prettier names
mapCytokineNamesFile = file.path('StimulusToPrettyNames.xlsx')
sheetName = 'cytokine'
loadedData = openxlsx::read.xlsx(xlsxFile = mapCytokineNamesFile, sheet = sheetName, colNames=F)
mappingListCyto = createMappingList(loadedData)

sheetName = 'circMed'
loadedData = openxlsx::read.xlsx(xlsxFile = mapCytokineNamesFile, sheet = sheetName, colNames=F)
mappingListCircMed = createMappingList(loadedData)

#Load peak annotations
peakAnnotations = fread(peakInfoFile, data.table = F)

#Define the ones we want to focus on
combinationsInMainPaper = list(
  c("cytoData.good.fc","V3"),
  c("cytoData.good.log2","V1"),
  c("circMarkerData.log2","V1"),
  c("cytoData.good.fc","V2")
)

combinationsInMainPaper.df = t(as.data.frame(combinationsInMainPaper))

cytokineGroups = defineSubgroupsCytokines()

createFigures= T
if (createFigures){
  ######Create the actual figures#########
  #============Create gene overview and make plots===========
  print(combinationsInMainPaper.df)
  for (rowIndex in c(1:nrow(combinationsInMainPaper.df))){
    #Define data and timepoint
    dataTypeSel=combinationsInMainPaper.df[[rowIndex,1]]
    timepointSel=combinationsInMainPaper.df[[rowIndex,2]]
    
    dataSubsetNamePretty = paste(dataTypeNameMapSupplement[[dataTypeSel]],mapVisitShortNameToFinal[[timepointSel]],sep='_')
    if (endsWith(suffix = '.fc', dataTypeSel)){
      dataSubsetNamePretty = paste0(dataSubsetNamePretty,'FC')
    }
    
    #Select data
    dataSel = allRelevantDataList[[dataTypeSel]]
    selVariables = colnames(dataSel)
    
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
    
    combinedResults <- loadAllAtacSeqAssociationResults(selFolder)

    subsetSel = 'GENE_AND_DISTAL_10kb'
    #Load the gene set data
    fileWithRelevantPeaks = file.path(folderPathways,paste0('peaks_filtered_PBMC.',subsetSel,'.csv.gz'))
    peakFileSel = fread(fileWithRelevantPeaks, data.table = F)
    
    combinedResultsSel = combinedResults[which(combinedResults$V1 %in% peakFileSel$peak_id),]
    
    #Ajust the pvalues
    combinedResultsSel$adj.P.Val.total = p.adjust(combinedResultsSel$p.value, method = 'fdr')
    combinedResultsSel$logPVal = -log10(combinedResultsSel$p.value)
    
    combinedResultsSel$'namesToPlot' = combinedResultsSel$gene_name
    combinedResultsSel$'namesToPlotItalic' = paste0("italic('",combinedResultsSel$gene_name,"')")
    combinedResultsSel = combinedResultsSel[order(combinedResultsSel$p.value, decreasing = F),]
    
    #For the fold-changes use more lenient cutoffs than for the absolue values
    if (endsWith(x = dataTypeSel, suffix = '.fc')){
      cutOffsSel = c(.25)
    }else{
      cutOffsSel = c(.1,.2,.3)
    }
    
    cutOffsSel = cutOffsSel[which(cutOffsSel>= min(combinedResultsSel$adj.P.Val.total))]
    
    if (length(cutOffsSel)==0){next}
    
    textSize = 12

    #Make the plots. Annotate the top hits in two ways: (1) With the top X, (2) top X plus top 1 per cytokines c(T,F)
    for (annotateTopPerCyto in c(F)){
      if (selData == 'cytokines'){
        plotOutput = plotBoxesWithTopHits(resultDf = combinedResultsSel,fdrCutOffs = cutOffsSel, minNrToAnnotate = 15, annotateTopPerCyto = annotateTopPerCyto, textSize=textSize)
      }else{
        plotOutput = plotBoxesWithTopHits(resultDf = combinedResultsSel,fdrCutOffs = cutOffsSel, minNrToAnnotate = 15, color.by = 'other', 
                                          annotateTopPerCyto = annotateTopPerCyto, textSize=textSize)
      }
      p = plotOutput[['p']]
      addToName = plotOutput[['addToName']]
      
      folderToSavePlot = file.path(image_dir_300BCG,'manuscriptImagesRaw','atacSeq','boxplots',paste(dataSubsetNamePretty,sep='_'))
      mkdirs(folderToSavePlot)
      
      sizeCombis=list()
      sizeCombis[["1"]] = c(2.5, 1.5, 6, 4.47)
      sizeCombis[["2"]] = c(6, 4.47, 6, 4.47)
      sizeCombis[["3"]] = c(6, 2*5.1216*0.393701, 6, 4.47)
      
      for (indexSize in names(sizeCombis)){
        print(indexSize)
        widthSel= sizeCombis[[indexSize]][[1]]#2
        widthPanel = sizeCombis[[indexSize]][[2]]#1.2
        heightSel = sizeCombis[[indexSize]][[3]]#5
        heightPanel = sizeCombis[[indexSize]][[4]]#4.47
        
        library(egg)
        #Fix the dimensions of the panel, so all of them nicely align
        p_fixed <- set_panel_size(p,
                                  width  = unit(widthPanel, "in"),
                                  height = unit(heightPanel, "in"))
        
        CairoFonts(regular="Helvetica:style=Regular",
                   bold="Helvetica:style=Bold",
                   italic="Helvetica:style=Italic",
                   bolditalic="Helvetica:style=Bold Italic,BoldItalic")
        ggsave(p_fixed, filename = file.path(folderToSavePlot,
                                             paste0(dataSubsetNamePretty,'_',subsetSel,'_',addToName,'_widthHeight=',widthSel,';',heightSel,'_',indexSize,'_new.pdf')),
               width=widthSel,height=heightSel, units = "in", device = cairo_pdf)
      }
    }
  }
}