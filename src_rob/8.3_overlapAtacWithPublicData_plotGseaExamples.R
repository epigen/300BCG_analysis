library(jsonlite)
library(R.utils)
library(data.table)
library(stringr)
library(dplyr)
library(openxlsx)
library(ggplot2)
source("plottingFxnsGeneral.R")

plotEnrichmentRaster <- function(pathway, stats, gseaParam = 1, ticksSize = 0.05, alphaSel=.2){
  library(fgsea)
  library(ggrastr)
  rnk <- rank(-stats)
  ord <- order(rnk)
  statsAdj <- stats[ord]
  statsAdj <- sign(statsAdj) * (abs(statsAdj)^gseaParam)
  statsAdj <- statsAdj/max(abs(statsAdj))
  pathway <- unname(as.vector(na.omit(match(pathway, names(statsAdj)))))
  pathway <- sort(pathway)
  gseaRes <- calcGseaStat(statsAdj, selectedStats = pathway, 
                          returnAllExtremes = TRUE)
  bottoms <- gseaRes$bottoms
  tops <- gseaRes$tops
  n <- length(statsAdj)
  xs <- as.vector(rbind(pathway - 1, pathway))
  ys <- as.vector(rbind(bottoms, tops))
  toPlot <- data.frame(x = c(0, xs, n + 1), y = c(0, ys, 0))
  diff <- (max(tops) - min(bottoms))/8
  x = y = NULL
  dpiSel=600
  p1 <- ggplot(toPlot, aes(x = x, y = y)) + 
    geom_hline(yintercept = max(tops), colour = "red", linetype = "dashed") + 
    geom_hline(yintercept = min(bottoms), colour = "red", linetype = "dashed") + 
    geom_hline(yintercept = 0, colour = "black") + 
    rasterise(geom_point(color = "green", size = 0.05), dpi=dpiSel) + 
    rasterise(geom_line(color = "green"), dpi=dpiSel) + 
    rasterise(geom_segment(data = data.frame(x = pathway), 
                           mapping = aes(x = x, y = -diff/2, xend = x, yend = diff/2), 
                           size = ticksSize, alpha = alphaSel, color='#000000'), dpi=dpiSel) + 
    theme_bw(base_size = 6, base_family = 'Helvetica') + 
    theme(panel.border = element_blank(), panel.grid.minor = element_blank()) + 
    theme(text=element_text(size=6,family = 'Helvetica')) +
    labs(x = "rank", y = "enrichment score") 
  p1

  return(p1)
}

###Get the variables we need
library(config)
Sys.setenv(R_CONFIG_ACTIVE = "300bcg")
config <- config::get()
for (x in names(config)){eval(parse(text=paste0(x,"='",config[[x]],"'")))}

folderToSaveImages = file.path(data_dir_300BCG,"atacSeq","trainedImmGeneSets_overlapResults","figures")
mkdirs(folderToSaveImages)

dirWithPublicResults = file.path(data_dir_300BCG, "atacSeq", "trainedImmGeneSets")
publicMetadata = fread(file = file.path(dirWithPublicResults,'processed','metadata_trainedImmunityGenesRegions.csv'),data.table = F)
publicMetadata$'public_data_folder' = file.path(dirWithPublicResults,'processed','converted_r97')

###List all the files that have results
basicSavePath = file.path(data_dir_300BCG,"atacSeq","trainedImmGeneSets_overlapResults")
allRankFiles = list.files(path = basicSavePath, pattern = 'ranks.*\\.RDS')
filePatterns = gsub('^ranks\\_\\_(.*)\\.RDS$','\\1',allRankFiles)

filePatternsSel = filePatterns

filePatternsSel = c("d90FC_Resp__main_2022.02.21.480081v1_5__top1000__up__distal",
                    "d90FC_Resp__main_29328908_6__top1000__up__distal",
                    "d0_NonR_vs_Resp__main_2022.02.21.480081v1_2__top1000__down__distal",
                    "d14FC_Resp__main_2022.02.21.480081v1_3__top1000__down__distal")

for (outputNameSel in filePatternsSel){
  
  if (grepl('.*\\_\\_main\\_.',outputNameSel)){
    metadataNameOtherData = 'publicMetadata'
  }else{
    metadataNameOtherData = 'scMetadata'
  }
  
  fgseaRes = readRDS(file.path(basicSavePath,
                               paste0(paste('fgseaRes',outputNameSel,sep='__'),'.RDS')))
  ranks = readRDS(file.path(basicSavePath,
                            paste0(paste('ranks',outputNameSel,sep='__'),'.RDS')))
  pathwaysListSingleEntry = readRDS(file.path(basicSavePath,
                                              paste0(paste('pathwaysListSingleEntry',outputNameSel,sep='__'),'.RDS')))
  
  comparisonSel = str_split(outputNameSel, pattern = '__', n = 2)[[1]][[1]]

  p1 = plotEnrichmentRaster(pathwaysListSingleEntry[[1]], ranks,
                            gseaParam = 1, ticksSize = 0.12,alphaSel=.1) #+ ggtitle(paste0('signed pval = ', pValSel))
  folderToSaveImage = file.path(folderToSaveImages,metadataNameOtherData,'GSEA')
  
  mkdirs(folderToSaveImage)
  ggsave(p1, 
         filename = file.path(folderToSaveImage,paste0(paste('GSEA',outputNameSel,sep='__'),'.pdf')),
         width=2.25*(3/4),height=1.4*(3/4))
}