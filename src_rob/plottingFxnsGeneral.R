
set_panel_size <- function(p=NULL, g=ggplotGrob(p), width=unit(3, "cm"), height=unit(3, "cm")){
  require(grid)
  require(ggplot2)
  panel_index_w<- g$layout$l[g$layout$name=="panel"]
  panel_index_h<- g$layout$t[g$layout$name=="panel"]
  g$widths[[panel_index_w]] <- width
  g$heights[[panel_index_h]] <- height
  class(g) <- c("fixed", class(g), "ggplot")
  g
}

print.fixed <- function(x) grid.draw(x)

left_width <- function(g){
  axis_l_index <- g$layout$r[g$layout$name=="axis-l"]
  ylab_index <- g$layout$r[g$layout$name=="ylab"]
  g$widths[[axis_l_index]] + g$widths[[ylab_index]]
}

full_width <- function(g){
  sum(g$widths)
}

align_plots <- function(..., width=unit(3, "cm"), height=unit(1, "null")){
  
  pl <- list(...)
  gl <- lapply(pl, set_panel_size, width=width, height=height)
  
  left <- lapply(gl, left_width)
  max_left <- max(do.call(unit.c, left))
  
  widths <- lapply(gl, full_width)
  max_width <- max(do.call(unit.c, widths))
  
  lay <- grid.layout(nrow=length(gl), ncol=1)
  vp <- viewport(layout=lay)
  pushViewport(vp)
  
  for(ii in seq_along(gl)){
    pushViewport(viewport(layout.pos.row=ii))
    pushViewport(viewport(x=unit(0.5, "npc") - 0.5*max_width + max_left - left[[ii]],
                          just="left", width=widths[[ii]]))
    grid.draw(gl[[ii]])
    upViewport(2)
  }
  upViewport()
}


sortMeasurementsForPlots <- function(measNames){
  subStrings = c('Macrophage_24h','PBMC_24h','PBMC_7d','Whole-Blood_48h','Ig')
  measNamesSort = c()
  for (subString in subStrings){
    indicesBool = grepl(subString, measNames)
    measNamesSortPart = sort(measNames[indicesBool])
    print(measNamesSortPart)
    measNamesSort = c(measNamesSort,measNamesSortPart)
  }
  measNamesSort = c(measNamesSort,sort(setdiff(measNames,measNamesSort)))
  return(measNamesSort)
}

plotHeatmapSignificance <- function(mm, sortedYAxis=F, sortedXAxis=F, image_dir, regressionType = 'regMethod', titleSel = 'title', 
                                    chosenRatio=1, scalingFactor=1, mm_months = F,looseFdrValues=F,plotFlipped=F,subdirsave='',
                                    extraColors = c(),extraColorLabels=c(), extraColorValues= c(), xAngle=30, use_cairo=F){
  #Rob ter Horst - Radboudumc - robterhorst88@gmail.com
  #mm - input data matrix
  #sortedYAxis - order to plot Y axis, if F will be done alphabetical
  #sortedXAxis - order to plot X axis, if F will be done alphabetical
  #image_dir - location to save
  #regressionType - type of regression used, used to add to name of file
  #titleSel - title
  #chosenRatio - ratio of the rectangles
  #mm_months - for plotting month names in the heamtap
  #looseFdrValues -  by default the highest signficant value is 0.05, if this parameters is set T the new higest colored pvalue is .1
  #plotFlipped - flip the heatmap (i.e. change X and Y)
  #subdirsave - 
  #extraColors - 2018-11-23 adding this ppart to allow extension of the palette with extra colors and values to have some extra indications
  #extraColorLabels - 2018-11-23 adding this ppart to allow extension of the palette with extra colors and values to have some extra indications
  #extraColorValues - the values in the df used for these values. The values should always be 10 or larger
  #extraColors = c("D08504","003031"),extraColorLabels=c('bothPos','bothNeg'), extraColorValues = c(10,20)
  
  library(scales)
  #library(gplots)
  library(RColorBrewer)
  library(ggplot2)
  
  if (use_cairo){library(Cairo)}

  stopifnot(all(extraColorValues>=10))
  
  if (plotFlipped){
    mm = t(mm)
    sortedYAxisTemp=sortedYAxis
    sortedYAxis=sortedXAxis
    sortedXAxis=sortedYAxisTemp
    regressionType = paste(regressionType,'flipped',sep='_')
  }
  is.nan.data.frame <- function(x){do.call(cbind, lapply(x, is.nan))}
  

  #mypalette <- c("#CECECE","#9b9b9b",brewer.pal(9,"Blues"),brewer.pal(9,"Reds"))
  mypalette <- c("#CECECE","#FFFFFF",brewer.pal(9,"Blues")[3:9],brewer.pal(9,"Reds")[3:9],extraColors)
  
  #Save the original data to use for the extra labels later
  mmRaw = mm
  
  if (looseFdrValues==1){
    mm[which(abs(mm)>0.1)]=NA
  }else if (looseFdrValues==2){
    mm[which(abs(mm)>0.2)]=NA
  }else{
    mm[which(abs(mm)>0.05)]=NA
  }
  
  mm.classes = mm
  mm.classes[abs(mm)>10] = NA
  
  if (looseFdrValues==1){
    mm.classes[mm>=-0.1 & mm<0] = 2
    mm.classes[mm>=-0.05 & mm<0] = 3
    mm.classes[mm>=-0.01 & mm<0] = 4
    mm.classes[mm>=-0.005 & mm<0] = 5
    mm.classes[mm>=-0.001 & mm<0] = 6
    mm.classes[mm>=-0.0001 & mm<0] = 7
    mm.classes[mm>=-0.00001 & mm<0] = 8
    mm.classes[mm<=0.1 & mm>0] = 9
    mm.classes[mm<=0.05 & mm>0] = 10
    mm.classes[mm<=0.001 & mm>0] = 11
    mm.classes[mm<=0.005 & mm>0] = 12
    mm.classes[mm<=0.001 & mm>0] = 13
    mm.classes[mm<=0.0001 & mm>0] = 14
    mm.classes[mm<=0.00001 & mm>0] = 15
    pValLabels = c("pre-filtered/not meas","not signif","<0.1 (neg)","<0.05 (neg)","<0.01 (neg)","<0.005 (neg)",
                   "<0.001 (neg)","<0.0001 (neg)","<0.00001 (neg)",
                   "<0.1 (pos)","<0.05 (pos)","<0.01 (pos)","<0.005 (pos)",
                   "<0.001 (pos)","<0.0001 (pos)","<0.00001 (pos)")
    
  }else if (looseFdrValues==2){
    mm.classes[mm>=-0.2 & mm<0] = 2
    mm.classes[mm>=-0.1 & mm<0] = 3
    mm.classes[mm>=-0.05 & mm<0] = 4
    mm.classes[mm>=-0.01 & mm<0] = 5
    mm.classes[mm>=-0.005 & mm<0] = 6
    mm.classes[mm>=-0.001 & mm<0] = 7
    mm.classes[mm>=-0.0001 & mm<0] = 8
    mm.classes[mm<=0.2 & mm>0] = 9
    mm.classes[mm<=0.1 & mm>0] = 10
    mm.classes[mm<=0.05 & mm>0] = 11
    mm.classes[mm<=0.01 & mm>0] = 12
    mm.classes[mm<=0.005 & mm>0] = 13
    mm.classes[mm<=0.001 & mm>0] = 14
    mm.classes[mm<=0.0001 & mm>0] = 15
    pValLabels = c("pre-filtered/not meas","not signif","<0.2 (neg)","<0.1 (neg)","<0.05 (neg)","<0.01 (neg)",
                   "<0.005 (neg)","<0.001 (neg)","<0.0001 (neg)",
                   "<0.2 (pos)","<0.1 (pos)","<0.05 (pos)","<0.01 (pos)",
                   "<0.005 (pos)","<0.001 (pos)","<0.0001 (pos)")
  }else{
    mm.classes[mm>=-0.05 & mm<0] = 2
    mm.classes[mm>=-0.01 & mm<0] = 3
    mm.classes[mm>=-0.005 & mm<0] = 4
    mm.classes[mm>=-0.001 & mm<0] = 5
    mm.classes[mm>=-0.0001 & mm<0] = 6
    mm.classes[mm>=-0.00001 & mm<0] = 7
    mm.classes[mm>=-0.000001 & mm<0] = 8
    mm.classes[mm<=0.05 & mm>0] = 9
    mm.classes[mm<=0.01 & mm>0] = 10
    mm.classes[mm<=0.005 & mm>0] = 11
    mm.classes[mm<=0.001 & mm>0] = 12
    mm.classes[mm<=0.0001 & mm>0] = 13
    mm.classes[mm<=0.00001 & mm>0] = 14
    mm.classes[mm<=0.000001 & mm>0] = 15
    pValLabels = c("pre-filtered/not meas","not signif","<0.05 (neg)","<0.01 (neg)","<0.005 (neg)","<0.001 (neg)",
                   "<0.0001 (neg)","<0.00001 (neg)","<0.000001 (neg)",
                   "<0.05 (pos)","<0.01 (pos)","<0.005 (pos)","<0.001 (pos)",
                   "<0.0001 (pos)","<0.00001 (pos)","<0.000001 (pos)")
  }
  
  pValLabels = c(pValLabels,extraColorLabels)
  if (length(extraColors)>0){
    labelnr=16
    for (extraColorIndex in 1:length(extraColorValues)){
      extraColorVal= extraColorValues[[extraColorIndex]]
      mm.classes[mmRaw==extraColorVal] = labelnr
      labelnr = labelnr + 1
    }
  }
  
  mm.classes[is.na(mm) & mmRaw<=1] = 1
  mm.classes[is.nan(mm)] = 0
  
  #color scale
  mm.classes.m = melt(as.matrix( mm.classes))
  colnames(mm.classes.m) = c('factorInfl','stimulus_cyto','value')
  #the ggplot
  base_size <- 15
  
  indicesToKeep = c(0:30) %in% unique(mm.classes.m$value)
  pValLabels = pValLabels[indicesToKeep]
  mypalette = mypalette[indicesToKeep]
  
  if (sortedYAxis==F){
    sortedYAxis <- sortMeasurementsForPlots(colnames(mm))
  }
  if (sortedXAxis== F){
    sortedXAxis <- sort(rownames(mm))
  }
  
  if (xAngle==90){
    hjust = 1
    vjust = .5
  }else{
    hjust = 1
    vjust = 1
  }
  
  #pValLabels= pValLabels[1:(max(mm.classes)+1)]
  p<-ggplot(mm.classes.m,aes(x=factorInfl,y=stimulus_cyto))+
    #tile layer
    geom_tile(aes(fill=factor(value)), colour = "black", show.legend=T, size=.75) + #, drop=F) +
    #setting the color
    scale_fill_manual(values=mypalette, 
                      labels=pValLabels, drop=FALSE) + 
    theme_grey(base_size = base_size) + labs(x = "", y = "") + scale_x_discrete(expand = c(0, 0)) + 
    theme(axis.ticks = element_blank(), 
          axis.text.x = element_text(size = 22, angle = xAngle, hjust = hjust, vjust = vjust, colour = "black"),
          axis.text.y = element_text(size = 22, colour = "black"),
          legend.text = element_text(size=22),
          legend.title = element_text(size=22),
          plot.title = element_text(hjust = 0.5, size = 22, colour = "black", face='bold'),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          panel.background = element_blank()
    ) +
    guides(fill=guide_legend(title="p-values")) +
    scale_y_discrete(limits=sort(colnames(mm.classes), decreasing =  TRUE), expand = c(0, 0)) +
    ylim(rev(sortedYAxis)) +  xlim(sortedXAxis) +
    ggtitle(titleSel)
  
  if (class(mm_months)=='data.frame'){
    mm_months.m = melt(as.matrix( mm_months))
    colnames(mm_months.m) = c('stimulus_cyto','factorInfl','value')
    
    theta <- seq(pi/8, 2*pi, length.out=64)
    xo <- .025
    yo <- .08
    for(i in theta) {
      p <- p + geom_text(data=mm_months.m, aes(x=factorInfl,y=stimulus_cyto,label=value), size=7, colour='black', hjust=.5+(cos(i)*xo), vjust =.5+(sin(i)*yo))
    }
    
    p= p +  geom_text(data=mm_months.m, size=7, colour = 'white', aes(x=factorInfl,y=stimulus_cyto, label = value)) 
  }
  
  chosenHeightPerBlock = 1.4 * scalingFactor
  chosenHeight = chosenHeightPerBlock*length(sortedYAxis)
  chosenWidth = chosenHeightPerBlock*length(sortedXAxis)*(1/chosenRatio)
  p=set_panel_size(p, width=unit(chosenWidth, "cm"), height = unit(chosenHeight, "cm")) #+ coord_fixed(ratio=chosenRatio)
  
  saveFolder = file.path(image_dir,'heatmaps',subdirsave)
  mkdirs(saveFolder)
  fileName = file.path(saveFolder,paste('heatmapSignificantFactors',Sys.time(),'_',regressionType,'.pdf',sep=''))
  plotWidth = chosenWidth*4 + (max(unlist(lapply(colnames(mm), nchar)))/2)
  plotHeight = chosenHeight*3 + 25#*2+7
  fileName = gsub('(\\-|// |\\:)','_',fileName)
  print(fileName)
  if (use_cairo){
    ggsave(plot = p, filename =  fileName, width = plotWidth, height=plotHeight, units = c("cm"), dpi = 300,limitsize=FALSE, device=cairo_pdf)
  }else{
    ggsave(plot = p, filename =  fileName, width = plotWidth, height=plotHeight, units = c("cm"), dpi = 300,limitsize=FALSE)#, device=cairo_pdf)
  }
  return(p)
}



plotCorrelationsHeatmap <- function(correlations, sortedYAxis=F, sortedXAxis=F, image_dir, addToName = 'regMethod', titleSel = 'title', 
                                    chosenRatio=1, scalingFactor=1, plotFlipped=F,subdirsave='', xAngle=30){
  #Rob ter Horst - Radboudumc - robterhorst88@gmail.com
  #correlations - input data matrix
  #sortedYAxis - order to plot Y axis, if F will be done alphabetical
  #sortedXAxis - order to plot X axis, if F will be done alphabetical
  #image_dir - location to save
  #addToName - type of regression used, used to add to name of file
  #titleSel - title
  #chosenRatio - ratio of the rectangles
  #plotFlipped - flip the heatmap (i.e. change X and Y)
  #subdirsave - 
  
  library(scales)
  library(RColorBrewer)
  library(Cairo)
  library(ggplot2)
  
  if (plotFlipped){
    correlations = t(correlations)
    sortedYAxisTemp=sortedYAxis
    sortedYAxis=sortedXAxis
    sortedXAxis=sortedYAxisTemp
    addToName = paste(addToName,'flipped',sep='_')
  }
  is.nan.data.frame <- function(x){do.call(cbind, lapply(x, is.nan))}
  
  if (sortedYAxis==F){
    sortedYAxis <- sortMeasurementsForPlots(colnames(correlations))
  }
  if (sortedXAxis== F){
    sortedXAxis <- sort(rownames(correlations))
  }
  
  correlations.m = melt(as.matrix( correlations))
  colnames(correlations.m) = c('factorInfl','stimulus_cyto','value')
  
  base_size <- 15

  p<-ggplot(correlations.m,aes(x=factorInfl,y=stimulus_cyto))+
    #tile layer
    geom_tile(aes(fill=value), colour = 'black',size=.75, show.legend=T) + 
    scale_fill_gradientn(colours=c(brewer.pal(9,"Blues")[[9]],"#FFFFFF",brewer.pal(9,"Reds")[[9]]),
                         guide = guide_colorbar(label = TRUE,
                                                draw.ulim = TRUE, 
                                                draw.llim = TRUE,
                                                # here comes the code change:
                                                frame.linewidth = 2,
                                                frame.colour = "black",
                                                ticks = TRUE, 
                                                barwidth = 1.3,
                                                barheight = 6.5, 
                                                direction = 'vertical')) + 
    theme_grey(base_size = base_size) + labs(x = "", y = "") + scale_x_discrete(expand = c(0, 0)) + 
    theme(axis.ticks = element_blank(), 
          axis.text.x = element_text(size = 22, angle = xAngle, hjust = 1, vjust = 1, colour = "black"),
          axis.text.y = element_text(size = 22, colour = "black"),
          legend.text = element_text(size=22),
          legend.title = element_text(size=22),
          plot.title = element_text(hjust = 0.5, size = 22, colour = "black", face='bold'),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          panel.background = element_blank(),
          legend.background = element_rect(colour = 'black', fill = 'white', linetype='solid')
    ) +
    ylim(rev(sortedYAxis)) +  xlim(sortedXAxis) +
    ggtitle(titleSel)
  
  
  chosenHeightPerBlock = 1.4 * scalingFactor
  chosenHeight = chosenHeightPerBlock*length(sortedYAxis)
  chosenWidth = chosenHeightPerBlock*length(sortedXAxis)*(1/chosenRatio)
  p=set_panel_size(p, width=unit(chosenWidth, "cm"), height = unit(chosenHeight, "cm"))
  
  saveFolder = file.path(image_dir,'heatmaps',subdirsave)
  mkdirs(saveFolder)
  fileName = file.path(saveFolder,paste('correlations',format(Sys.time(), format = "%Y-%m-%d_%H-%M-%S"),'_',addToName,'.pdf',sep=''))
  plotWidth = chosenWidth*4 + (max(unlist(lapply(colnames(correlations), nchar)))/2)
  plotHeight = chosenHeight*3 + 25#*2+7
  print(fileName)
  ggsave(plot = p, filename =  fileName, width = plotWidth, height=plotHeight, units = c("cm"), dpi = 300,limitsize=FALSE, device=cairo_pdf)
}