library(jsonlite)
library(R.utils)
library(data.table)
library(stringr)
library(dplyr)
library(openxlsx)
library(ggplot2)
source("plottingFxnsGeneral.R")


###Get the variables we need
library(config)
Sys.setenv(R_CONFIG_ACTIVE = "300bcg")
config <- config::get()
for (x in names(config)){eval(parse(text=paste0(x,"='",config[[x]],"'")))}

fgseaPerm=10^6

folderToSaveOverlapResults = file.path(data_dir_300BCG,"atacSeq","trainedImmGeneSets_overlapResults")

###ATAC results 
dirWithAtacResults = file.path(data_dir_300BCG, "atacSeq","trainedImmGeneSets","300BCG_results")
resultFiles = list.files(dirWithAtacResults,pattern = '^TRIM.*')
comparisonsToCheck = gsub('^TRIM\\.(.*)\\.csv$','\\1',resultFiles)

#Other ATAC results 
resultFiles2 = list.files(dirWithAtacResults,pattern = '^PC2.*')
comparisonsToCheck2 = gsub('\\.csv','',resultFiles2)

resultFilesFull = c(resultFiles,resultFiles2)
comparisonsToCheckFull = c(comparisonsToCheck,comparisonsToCheck2)

###Processed public data
dirWithPublicResults = file.path(data_dir_300BCG, "atacSeq", "trainedImmGeneSets")
publicMetadata = fread(file = file.path(dirWithPublicResults,'processed','metadata_trainedImmunityGenesRegions.csv'),data.table = F)
publicMetadata$'public_data_folder' = file.path(dirWithPublicResults,'processed','converted_r97')

addScData=F
if (addScData){
  ###Now add the single cell results in the same framework as the public data
  singleCellFiles = list.files(dirWithAtacResults,pattern = '^scRNAseq.scLM.*')
  nrowVec = c()
  for (fileSel in singleCellFiles){
    scData = fread(file = file.path(dirWithAtacResults,singleCellFiles[[1]]))
    nrowVec = c(nrowVec,nrow(scData))
  }
  
  singleCellComparisons = gsub('^scRNAseq\\.(.*)\\.csv','\\1',singleCellFiles)
  scMetadata = as.data.frame(matrix(NA,ncol=length(publicMetadata),nrow = length(singleCellComparisons)))
  colnames(scMetadata) = colnames(publicMetadata)
  scMetadata$filename = gsub('\\.csv$','',singleCellFiles)
  scMetadata$assay = "scRNA-seq" 
  scMetadata$direction = "both" 
  scMetadata$PMID = "privateScData"
  scMetadata$comparison = singleCellComparisons
  scMetadata$cell_type = "PBMC"
  scMetadata$orig_sheet = "not applicable"
  scMetadata$full_result = "yes"
  scMetadata$description = singleCellComparisons
  scMetadata$orig_columns = "human_gene_ensembl"
  scMetadata$data_type = "gene scores"
  scMetadata$organism = "hs"
  scMetadata$region_filter = "enhancers"
  scMetadata[grep('TSS_PROXIMAL',singleCellComparisons),'region_filter'] = "promoters"
  scMetadata$orig_filename = singleCellFiles
  scMetadata$n_rows = nrowVec
  scMetadata$'public_data_folder' = dirWithAtacResults
}

folderToSaveImages = file.path(data_dir_300BCG,"atacSeq","trainedImmGeneSets_overlapResults","figures")
mkdirs(folderToSaveImages)

if (addScData){
  metadataList = list('publicMetadata'=publicMetadata, 'scMetadata'=scMetadata)
}else{
  metadataList = list('publicMetadata'=publicMetadata)
}

comparisonsToCheckSel = comparisonsToCheckFull
comparisonsToCheckSel[[which(comparisonsToCheckSel=="d0_Resp_vs_NonR")]] = "d0_NonR_vs_Resp" #Convention we changed for visualization in manuscript
for (metadataNameOtherData in names(metadataList)){
  metadataOtherData = metadataList[[metadataNameOtherData]]
  ##Now combine the data and make the plots in the combinations we want
  studiesToCheck = metadataOtherData$filename
  resultsDfPerCombiEmpty = as.data.frame(matrix(NA,ncol = length(comparisonsToCheckSel), nrow=length(studiesToCheck)))
  colnames(resultsDfPerCombiEmpty) = comparisonsToCheckSel
  rownames(resultsDfPerCombiEmpty) = studiesToCheck
  
  basicComparison = 'allReported'
  regionUniverses = c('promoter','distal')
  
  resultsDfList = list()
  for (testTypeSel in c('fisher','GSEA')){#
    for (topSel in c('top1000','top2000','p_0_05','padj_0_1')){
      for (regionUniverseSel in c('promoter','distal')){
        for (directionSel in c('up','down')){
          
          resultsDfCombi = resultsDfPerCombiEmpty
          for (comparisonSel in rev(comparisonsToCheckSel)){ #c("d90FC_Resp","d14FC_Resp")){
            
            
            for (filenameSel in metadataOtherData[,'filename']){
              
              fileToLoad = file.path(data_dir_300BCG,"atacSeq","trainedImmGeneSets_overlapResults",
                                     paste0(comparisonSel,'__',metadataNameOtherData,"__",filenameSel,'__nPermFgsea=',fgseaPerm,".RDS"))
              if (!file.exists(fileToLoad)){next}
              
              resultsOutputList = readRDS(file = fileToLoad)
              if (regionUniverseSel %in% names(resultsOutputList)){
                if (directionSel %in% names(resultsOutputList[[regionUniverseSel]])){
                  
                  if (basicComparison %in% names(resultsOutputList[[regionUniverseSel]][[directionSel]])){
                    pValSel = resultsOutputList[[regionUniverseSel]][[directionSel]][[basicComparison]][[testTypeSel]]
                  }else{
                    pValSel = resultsOutputList[[regionUniverseSel]][[directionSel]][[topSel]][[testTypeSel]]
                  }
                  if (is.null(pValSel)){pValSel=NA}
                  
                  resultsDfCombi[[filenameSel,comparisonSel]] = pValSel
                }
              }
            }
          }
          
          resultsDfList[[paste(testTypeSel,regionUniverseSel,topSel,directionSel,sep='__')]] = resultsDfCombi
          
          
          matrixForPlot = as.matrix(resultsDfCombi)
          matrixForPlot[is.na(matrixForPlot)]=NaN
          
          folderToSaveHeatmaps = file.path(folderToSaveImages,metadataNameOtherData)
          mkdirs(folderToSaveHeatmaps)
          plotHeatmapSignificance(matrixForPlot, sortedYAxis=F, sortedXAxis=F, 
                                  image_dir = folderToSaveHeatmaps, 
                                  regressionType = paste(testTypeSel,regionUniverseSel,topSel,directionSel,sep='__'), 
                                  titleSel = paste(testTypeSel,regionUniverseSel,topSel,directionSel,sep='__'), 
                                  chosenRatio=1, scalingFactor=1, mm_months = F,looseFdrValues=F,plotFlipped=F,subdirsave='',
                                  extraColors = c(),extraColorLabels=c(), extraColorValues= c(), xAngle=90)
          
          ###Save matrix for supplement
          openxlsx::write.xlsx(as.data.frame(t(matrixForPlot)),
                               file.path(suppTableDir,'xlsx_small','raw',paste0('overlapAnalysis_',
                                                                                paste(metadataNameOtherData,testTypeSel,regionUniverseSel,topSel,directionSel,sep='__'),
                                                                                '_v1.xlsx')),
                               col.names=TRUE, row.names=TRUE)
        }
      }
    }
  }
  
  
  ###Now make some combined plots
  for (testTypeSel in c('GSEA')){
    for (topSel in c('top1000')){ #other options would be 'top2000','p_0_05','padj_0_1'
      for (regionUniverseSel in c('distal')){  #can also be set to 'promoter'
        resultsDf_up = resultsDfList[[paste(testTypeSel,regionUniverseSel,topSel,'up',sep='__')]]
        colnames(resultsDf_up) = paste(colnames(resultsDf_up),'up',sep='_')
        resultsDf_down = resultsDfList[[paste(testTypeSel,regionUniverseSel,topSel,'down',sep='__')]]
        colnames(resultsDf_down) = paste(colnames(resultsDf_down),'down',sep='_')
        stopifnot(all.equal(rownames(resultsDf_up), rownames(resultsDf_down)))
        
        resultsDf = cbind.data.frame(resultsDf_up, resultsDf_down)
        
        if (metadataNameOtherData=='scMetadata'){
          if (regionUniverseSel == 'distal'){
            resultsDf = resultsDf[grepl('GENE_AND_DISTAL_10kb',rownames(resultsDf)),]
          }else if (regionUniverseSel == 'promoter'){
            resultsDf = resultsDf[grepl('TSS_PROXIMAL',rownames(resultsDf)),]
          }
        }else if (metadataNameOtherData=='publicMetadata'){
          
          orderToPlotCellTypes = c("monocyte",
                                   "granulocyte-monocyte progenitor",
                                   "macrophage",
                                   "BMDM",
                                   'MDM',
                                   "alveolar-like MDM",
                                   "induced sputum-isolated macrophage",
                                   'BMDDC',
                                   "HSC",
                                   "HSPC",
                                   "HSC+MPP+GMP",
                                   "MPP",
                                   "MPP3",
                                   "MPP4",
                                   "various",
                                   "CD34+",
                                   "LSK"
          )
          
          #Now change rownames
          newRownames = rownames(resultsDf)
          for (rowIndex in c(1:nrow(resultsDf))){
            filename=rownames(resultsDf)[[rowIndex]]
            assaySel=publicMetadata[publicMetadata$filename==filename,'assay']
            organism = publicMetadata[publicMetadata$filename==filename,'organism']
            cell_type = publicMetadata[publicMetadata$filename==filename,'cell_type']
            
            ###Do some mappings of rownames
            #monocyte-derived macrophage = MDM
            #bone marrow–derived DC = BMDDC
            #aMDM = alveolar-like MDM
            #alveolar-like monocyte-derived macrophage = alveolar-like MDM
            
            mappingVector = c('MDM',
                              'BMDDC',
                              'alveolar-like MDM',
                              'alveolar-like MDM')
            names(mappingVector) = c('monocyte-derived macrophage',
                                     'bone marrow–derived DC',
                                     'aMDM',
                                     'alveolar-like monocyte-derived macrophage')
            if (cell_type %in% names(mappingVector)){
              cell_type_mapped = mappingVector[[cell_type]]
            }else{
              cell_type_mapped = cell_type
            }
            
            newRownames[[rowIndex]] = paste(
              assaySel,
              organism,
              cell_type_mapped,
              gsub('^main\\_','',rownames(resultsDf)[[rowIndex]]),
              sep= '  ')
            print(newRownames[[rowIndex]])
          }
          
          orderToPlotAssays = c(
            "ATAC-seq",
            "H3K4me1",
            "H3K4me3",
            "H3K27ac",
            "RNA-seq",
            "scRNA-seq",
            "Illumina Infinium MethylationEPIC BeadChips",
            "pathway_databases",
            "various"
          )
          
          orderToPlotOrganisms = c("hs","mm","various")
          
          rowInfo = as.data.frame(str_split_fixed(string = newRownames, pattern ='  ', n = 4))
          colnames(rowInfo) = c('assay','organism','cell_type','study')
          
          rowInfo$pmid = gsub(pattern = '\\_([0-9]+)$',replacement = '',x = as.character(rowInfo$study))
          rowInfo$nrInStudy = as.numeric(gsub(pattern = '.*_([0-9]+)$',replacement = '\\1',x = as.character(rowInfo$study)))
          
          rowInfo[['assay']] = factor(rowInfo[['assay']], levels=orderToPlotAssays)
          rowInfo[['cell_type']] = factor(rowInfo[['cell_type']], levels=orderToPlotCellTypes)
          rowInfo[['organism']] = factor(rowInfo[['organism']], levels=orderToPlotOrganisms)
          rowInfo$oldOrder = c(1:nrow(rowInfo))
          
          openxlsx::write.xlsx(rowInfo, file = file.path(folderToSaveImages,
                                                         metadataNameOtherData,
                                                         'rowInfo.xlsx'))
          
          newOrder = with(rowInfo, order(assay, organism, cell_type, pmid, nrInStudy))
          
          #Change rownames
          newRownames2 = paste(c(1:nrow(rowInfo)), rowInfo[,'assay'], rowInfo[,'cell_type'], rowInfo[,'organism'])
          rownames(resultsDf) = newRownames2
          resultsDf = resultsDf[newOrder,]
          rowInfo = rowInfo[newOrder,]
        }
        
        ##Change names to what is more standardized
        mappingComparisons = c("Non-R. vs. Resp.","Day 14 vs Day 0","Day 90 vs Day 0")
        names(mappingComparisons) = c("d0_NonR_vs_Resp","d14FC_Resp","d90FC_Resp")
        for (nameSel in names(mappingComparisons)){
          colnames(resultsDf) = gsub(nameSel,mappingComparisons[[nameSel]], colnames(resultsDf))
        }
        
        source("plottingFxnsGeneral.R")
        matrixForPlot = as.matrix(resultsDf)
        matrixForPlot[is.na(matrixForPlot)]=NaN
        directionSel='up_and_down_combined'
        
        folderToSaveHeatmaps = file.path(folderToSaveImages,metadataNameOtherData)
        mkdirs(folderToSaveHeatmaps)
        
        plotHeatmapSignificance(matrixForPlot, sortedYAxis=colnames(matrixForPlot), sortedXAxis=rownames(matrixForPlot), 
                                image_dir = folderToSaveHeatmaps, 
                                regressionType = paste(testTypeSel,regionUniverseSel,topSel,directionSel,sep='__'), 
                                titleSel = paste(testTypeSel,regionUniverseSel,topSel,directionSel,sep='__'), 
                                chosenRatio=1, scalingFactor=1, mm_months = F,looseFdrValues=F,plotFlipped=F,subdirsave='',
                                extraColors = c(),extraColorLabels=c(), extraColorValues= c(), xAngle=90)
        
        openxlsx::write.xlsx(as.data.frame(t(matrixForPlot)), 
                             file.path(suppTableDir,'xlsx_small','raw',paste0('overlapAnalysis_',
                                                                              paste(metadataNameOtherData,testTypeSel,regionUniverseSel,topSel,directionSel,sep='__'),
                                                                              '_v2.xlsx')),
                             col.names=TRUE, row.names=TRUE)
        
        ###Now also plot one where the results are clustered
        #First remove 
        matrixForPlot.filt = matrixForPlot[which(rowSums(is.na(matrixForPlot))<ncol(matrixForPlot)),]
        #Cluster
        matrixForPlot.filt.log10 = -log10(abs(matrixForPlot.filt)) * sign(matrixForPlot.filt)
        #
        correlations = cor(t(matrixForPlot.filt.log10), use='pairwise.complete.obs')
        correlations[is.na(correlations)] = 0
        clusteringRes = hclust(as.dist(1-correlations), method = "ward.D2")
        # Visualization using the default theme named theme_dendro()
        library(ggdendro)
        folderToSaveDendro= file.path(folderToSaveImages,metadataNameOtherData,"hclust")
        mkdirs(folderToSaveDendro)
        p2 = ggdendrogram(clusteringRes)
        ggsave(plot = p2, 
               filename = file.path(folderToSaveDendro,paste0(paste(testTypeSel,regionUniverseSel,topSel,directionSel,sep='__'),"_hclust.pdf")),
               width = 5, height = 6)
        
        matrixForPlot.filt = matrixForPlot.filt[clusteringRes$order,]
        matrixForPlot.filt = as.matrix(matrixForPlot.filt)
        matrixForPlot.filt[is.na(matrixForPlot.filt)]=NaN
        plotHeatmapSignificance(matrixForPlot.filt, sortedYAxis=colnames(matrixForPlot.filt), sortedXAxis=rownames(matrixForPlot.filt), 
                                image_dir = folderToSaveHeatmaps, 
                                regressionType = paste(testTypeSel,regionUniverseSel,topSel,directionSel,'clustered',sep='__'), 
                                titleSel = paste(testTypeSel,regionUniverseSel,topSel,directionSel,sep='__'), 
                                chosenRatio=1, scalingFactor=1, mm_months = F,looseFdrValues=F,plotFlipped=F,subdirsave='',
                                extraColors = c(),extraColorLabels=c(), extraColorValues= c(), xAngle=90)
        
        
        
        ##If PC2 is in there, make a separate plot for this one
        if (sum(grepl(pattern = '^PC2.*',colnames(matrixForPlot)))>0){
          ###Now also plot one where the results are clustered
          #First remove 
          matrixForPlot.filt = matrixForPlot[,which(grepl('^PC2.*',colnames(matrixForPlot)))]
          matrixForPlot.filt = matrixForPlot.filt[which(rowSums(is.na(matrixForPlot.filt))<ncol(matrixForPlot.filt)),]
          #
          print(head(matrixForPlot.filt))
          rowsToKeep=which(rowSums(abs(matrixForPlot.filt)<.01, na.rm=T)>=1)
          
          if (length(rowsToKeep)>2){
            
            matrixForPlot.filt = matrixForPlot.filt[rowsToKeep,]
            #Cluster
            matrixForPlot.filt.log10 = -log10(abs(matrixForPlot.filt)) * sign(matrixForPlot.filt)
            
            matrixForPlot.filt.log10.classes = matrixForPlot.filt.log10
            matrixForPlot.filt.log10.classes[]=0
            matrixForPlot.filt.log10.classes[abs(matrixForPlot.filt)<=.05 & sign(matrixForPlot.filt)>0]=1
            matrixForPlot.filt.log10.classes[abs(matrixForPlot.filt)<=.05 & sign(matrixForPlot.filt)<0]=-1
            
            matrixForPlot.filt.log10.classes2 = matrixForPlot.filt.log10.classes
            matrixForPlot.filt.log10.classes2[is.na(matrixForPlot.filt.log10.classes2)]=0
            distForClust = dist(matrixForPlot.filt.log10.classes2, method = "manhattan")
            #Add a small dist for the case where there is one NA
            distAsMatrix = as.matrix(distForClust)
            for (rowIndex in c(1:nrow(distAsMatrix))){
              for (colIndex in c(1:ncol(distAsMatrix))){
                if ((sum(is.na(matrixForPlot.filt[rowIndex,])) + sum(is.na(matrixForPlot.filt[colIndex,])))==1){
                  distAsMatrix[rowIndex,colIndex] = distAsMatrix[rowIndex,colIndex]+.5
                }
              }
            }
            
            clusteringRes = hclust(as.dist(distAsMatrix), method = "ward.D2")
            
            # Visualization using the default theme named theme_dendro()
            library(ggdendro)
            folderToSaveDendro= file.path(folderToSaveImages,metadataNameOtherData,"hclust")
            mkdirs(folderToSaveDendro)
            p2 = ggdendrogram(clusteringRes)
            ggsave(plot = p2, 
                   filename = file.path(folderToSaveDendro,paste0(paste(testTypeSel,regionUniverseSel,topSel,directionSel,sep='__'),"_hclust_PC2.pdf")),
                   width = 5, height = 6)
            
            matrixForPlot.filt = matrixForPlot.filt[clusteringRes$order,]
            
            matrixForPlot.filt = as.matrix(matrixForPlot.filt)
            matrixForPlot.filt[is.na(matrixForPlot.filt)]=NaN
            plotHeatmapSignificance(matrixForPlot.filt, sortedYAxis=colnames(matrixForPlot.filt), sortedXAxis=rownames(matrixForPlot.filt), 
                                    image_dir = folderToSaveHeatmaps, 
                                    regressionType = paste(testTypeSel,regionUniverseSel,topSel,directionSel,'clustered_PC2',sep='__'), 
                                    titleSel = paste(testTypeSel,regionUniverseSel,topSel,directionSel,sep='__'), 
                                    chosenRatio=1, scalingFactor=1, mm_months = F,looseFdrValues=F,plotFlipped=F,subdirsave='',
                                    extraColors = c(),extraColorLabels=c(), extraColorValues= c(), xAngle=90)
          }
        }
        
        ###Now split back up and keep only the significant when split into up -- down and epigenetics -- RNAseq
        if (metadataNameOtherData=='publicMetadata'){
          epigenAssays = c("ATAC-seq", "H3K4me1", "H3K4me3", "H3K27ac") 
          rnaAssays = c("RNA-seq", "scRNA-seq")
          epigenetics_up = resultsDf[which(rowInfo$assay %in% epigenAssays),c("Non-R. vs. Resp._up","Day 14 vs Day 0_up","Day 90 vs Day 0_up")]
          epigenetics_down = resultsDf[which(rowInfo$assay %in% epigenAssays),c("Non-R. vs. Resp._down","Day 14 vs Day 0_down","Day 90 vs Day 0_down")]
          rnaseq_up = resultsDf[which(rowInfo$assay %in% rnaAssays),c("Non-R. vs. Resp._up","Day 14 vs Day 0_up","Day 90 vs Day 0_up")]
          rnaseq_down = resultsDf[which(rowInfo$assay %in% rnaAssays),c("Non-R. vs. Resp._down","Day 14 vs Day 0_down","Day 90 vs Day 0_down")]
          splitDfList = list()
          splitDfList[["epigenetics_up"]] = epigenetics_up
          splitDfList[["epigenetics_down"]] = epigenetics_down
          splitDfList[["rnaseq_up"]] = rnaseq_up
          splitDfList[["rnaseq_down"]] = rnaseq_down
        }else{
          rnaseq_up = resultsDf[,c("Non-R. vs. Resp._up","Day 14 vs Day 0_up","Day 90 vs Day 0_up")]
          rnaseq_down = resultsDf[,c("Non-R. vs. Resp._down","Day 14 vs Day 0_down","Day 90 vs Day 0_down")]
          splitDfList = list()
          splitDfList[["rnaseq_up"]] = rnaseq_up
          splitDfList[["rnaseq_down"]] = rnaseq_down
        }
        
        ##Filter for just keeping ones where at least one is significant
        for (subgroupSel in names(splitDfList)){
          resultsDfSel = splitDfList[[subgroupSel]]
          #Remove the ones with only insignificant or missing
          rowsToKeep=which((rowSums(abs(resultsDfSel)>.01, na.rm=T) + rowSums(is.na(resultsDfSel)))<3)
          if (length(rowsToKeep)>=2){
            resultsDfSel = resultsDfSel[rowsToKeep,]
            resultsDfSel.log10 = -log10(abs(resultsDfSel)) * sign(resultsDfSel)
            
            clusteringRes = hclust(dist(resultsDfSel.log10), method = "ward.D2")
            
            # Visualization using the default theme named theme_dendro()
            library(ggdendro)
            folderToSaveDendro= file.path(folderToSaveImages,metadataNameOtherData,"hclust")
            mkdirs(folderToSaveDendro)
            p2 = ggdendrogram(clusteringRes)
            ggsave(plot = p2, 
                   filename = file.path(folderToSaveDendro,paste0(paste(testTypeSel,regionUniverseSel,topSel,subgroupSel,sep='__'),"_hclust.pdf")),
                   width = 5, height = 6)
            
            resultsDfSel = resultsDfSel[clusteringRes$order,]
            
            matrixForPlot = as.matrix(resultsDfSel)
            matrixForPlot[is.na(matrixForPlot)]=NaN
            plotHeatmapSignificance(matrixForPlot, sortedYAxis=colnames(matrixForPlot), sortedXAxis=rownames(matrixForPlot), 
                                    image_dir = folderToSaveHeatmaps, 
                                    regressionType = paste(testTypeSel,regionUniverseSel,topSel,subgroupSel,sep='__'), 
                                    titleSel = paste(testTypeSel,regionUniverseSel,topSel,subgroupSel,sep='__'), 
                                    chosenRatio=1, scalingFactor=1, mm_months = F,looseFdrValues=F,plotFlipped=F,subdirsave='',
                                    extraColors = c(),extraColorLabels=c(), extraColorValues= c(), xAngle=90)
          }
        }
      }
    }
  }
}