# x should be numeric
inverseRankNorm  <-  function(x){
  #Inverse rank based
  #Made by Raul raul.aguirre.gamboa@gmail.com
  res <- rank(x)
  res <- qnorm(res/(length(res)+0.5))
  return(res)
}

inverseRankNormDf <- function(dataFrameSel){
  #This does the inverse rank transform for a dataframe,
  #The function raul gave me does not allow for NA's,
  # i have to run the algorithm separately for just 
  #the data with no NA's and then add it back to the
  #full Df
  dataFrameSelPart = dataFrameSel
  dataFrameSelBu = dataFrameSel
  for (colNameSel in colnames(dataFrameSelPart)){
    selectedData = dataFrameSelPart[,colNameSel]
    notNaIndices = which(!is.na(selectedData))
    selectedData = selectedData[notNaIndices]
    selectedData =  inverseRankNorm(selectedData)
    dataFrameSel[notNaIndices,colNameSel] = selectedData
    #The ranks in the transformed data should be the same as in the original data
    stopifnot(all.equal(rank(dataFrameSel[,colNameSel]),rank(dataFrameSelBu[,colNameSel])))
  }
  return(dataFrameSel)
}

mapCytokineNames <- function(nameVector,dataType,src_dir, mappingFileName = 'GeneralMappingOfStimulusNames.xlsx'){
  print(paste0('mapping with ',mappingFileName))
  createMappingList <- function(loadedData){
    #Mapping list with values first column
    mappingList = list()
    for (rowIndex in 1:nrow(loadedData)){
      valueSel = loadedData[rowIndex,1]
      for (colSel in 2:ncol(loadedData)){
        keySel = loadedData[rowIndex,colSel]
        if (!is.na(keySel)){
          mappingList[[keySel]] = valueSel
        }
      }
    }
    return(mappingList)
  }
  
  mapCytokineNamesFile= file.path(mappingFileName)
  
  stopifnot(dataType %in% c('stimCyto','circMed','stim','cyto'))
  
  if (dataType == 'stimCyto'){
    #Load mapping info
    sheetName = 'stimulus'
    loadedData = openxlsx::read.xlsx(xlsxFile = mapCytokineNamesFile, sheet = sheetName, colNames=F)
    mappingListStim = createMappingList(loadedData)
    
    sheetName = 'cytokine'
    loadedData = openxlsx::read.xlsx(xlsxFile = mapCytokineNamesFile, sheet = sheetName, colNames=F)
    mappingListCyto = createMappingList(loadedData)
    
    #Now split the cyto data info and map
    cytoInfo = str_split_fixed(nameVector,'_',4)
    cytoInfoNew = cytoInfo
    loadedStims = cytoInfo[,1]
    loadedCytos = cytoInfo[,4]
    print(setdiff(loadedStims, names(mappingListStim)))
    print(setdiff(loadedCytos, names(mappingListCyto)))
    stopifnot(all(loadedStims %in% names(mappingListStim)))
    stopifnot(all(loadedCytos %in% names(mappingListCyto)))
    cytoInfoNew[,1] = unlist(lapply(loadedStims, function(x) mappingListStim[[x]]))
    cytoInfoNew[,4] = unlist(lapply(loadedCytos, function(x) mappingListCyto[[x]]))
    
    #New namevector
    newNameVector = apply(cytoInfoNew, 1, paste, collapse="_")
  } else if (dataType == 'circMed') {
    sheetName = 'circMed'
    loadedData = openxlsx::read.xlsx(xlsxFile = mapCytokineNamesFile, sheet = sheetName, colNames=F)
    mappingList = createMappingList(loadedData)
    stopifnot(all(nameVector %in% names(mappingList)))
    newNameVector = unlist(lapply(nameVector, function(x) mappingList[[x]]))
  } else if (dataType == 'stim') {
    sheetName = 'stimulus'
    loadedData = openxlsx::read.xlsx(xlsxFile = mapCytokineNamesFile, sheet = sheetName, colNames=F)
    mappingListStim = createMappingList(loadedData)
    stopifnot(all(nameVector %in% names(mappingListStim)))
    newNameVector = unlist(lapply(nameVector, function(x) mappingListStim[[x]]))
    
  }else if (dataType == 'cyto') {
    sheetName = 'cytokine'
    loadedData = openxlsx::read.xlsx(xlsxFile = mapCytokineNamesFile, sheet = sheetName, colNames=F)
    mappingListCyto = createMappingList(loadedData)
    stopifnot(all(nameVector %in% names(mappingListCyto)))
    newNameVector = unlist(lapply(nameVector, function(x) mappingListCyto[[x]]))
  }
  
  
  #print statement
  print(newNameVector)
  together= cbind(nameVector,newNameVector)
  print(together)
  return(newNameVector)
}

log10TransformWithChecks <- function(selDataSet){
  #First check for 0's
  stopBecauseOfZeros = F
  for (selVar in colnames(selDataSet)){
    if (sum(selDataSet[[selVar]]==0, na.rm=T)>0){
      stopBecauseOfZeros = T
      print(paste0('zeros in ',selVar,' so not suitable for log transform'))
    }
    
  }
  if (stopBecauseOfZeros){
    stop('there are zeros in the data you want to log-transform')
  }
  
  #Log transform
  selDataSetLog = log10(selDataSet)
  
  #Finally, check if any extra inf, -inf or NA NaN were indtroduced
  stopifnot(sum(is.infinite(as.matrix(selDataSetLog)))==sum(is.infinite(as.matrix(selDataSet))))
  
  return(selDataSetLog)
}

calculateResiduals <- function(matrixInput, covarsMatrix, selCovars){
  matrixResiduals = as.data.frame(matrix(NA, nrow = nrow(matrixInput),ncol = ncol(matrixInput)))
  colnames(matrixResiduals) = colnames(matrixInput)
  rownames(matrixResiduals) = rownames(matrixInput)
  for (dependentVar in colnames(matrixInput)){
    print(dependentVar)
    mA1 <- lm(data = cbind(matrixInput,covarsMatrix), formula = paste0(dependentVar, " ~ ", paste(selCovars, collapse = ' + ') ), na.action = na.exclude ) #Create a linear model
    matrixResiduals[,dependentVar] = resid(mA1) 
  }
  return(matrixResiduals)
}

completeDataToFullSetRetainClasses <- function(dataFrame, hvIds){
  #Create a data frame of he data in the correct order, with no missing hvIds,
  #missing values are supplemented with NAs
  #First check if there are any double identifiers
  countIdent = sort(table(rownames(dataFrame)))
  
  if (max(countIdent)>=2){
    stop(paste("some identifiers appear more than once:", paste(names(countIdent), collapse = ' '),sep=' '))
  }
 
  #Here I Add the NA's
  for (hvId in hvIds){
    if (!(hvId %in% rownames(dataFrame))){
      newRow = data.frame(matrix(NA, nrow = 1, ncol=dim(dataFrame)[2]))
      rownames(newRow)=hvId
      colnames(newRow)= colnames(dataFrame)
      dataFrame = rbind(dataFrame, newRow)
    }
  }
  #Here I order them
  dataFrameOrdered = dataFrame[hvIds,,drop=F]
  
  return(dataFrameOrdered)
}

calcExtremeThresholdsTuckey <- function(data, iqrThres){
  lowerq = quantile(data, na.rm = T)[2]
  upperq = quantile(data, na.rm = T)[4]
  iqr = upperq - lowerq #Or use IQR(data)
  extreme.threshold.upper = (iqr * iqrThres) + upperq
  extreme.threshold.lower = lowerq - (iqr * iqrThres)
  return(c(extreme.threshold.lower,extreme.threshold.upper))
}

#Tukey's methods for outlier detection
#http://stackoverflow.com/questions/12866189/calculating-the-outliers-in-r
removeOutliersTuckeys <- function(selDataSet, iqrThres=6){
  #print(paste('IQR threshold =', iqrThres))
  selDataSetOutliersRem = cbind(selDataSet)
  for (varName in colnames(selDataSet)){
    #print(varName)
    data = selDataSetOutliersRem[,varName]
    if (sum(is.na(data))>.75*length(data)){
      print(paste('too many NA"s for ',varName))
      next()
      Sys.sleep(2)
    }
    iqrThresholds = calcExtremeThresholdsTuckey(data, iqrThres)
    extreme.threshold.upper =  iqrThresholds[[2]]
    extreme.threshold.lower = iqrThresholds[[1]]
    if (sum(data<=extreme.threshold.lower, na.rm = T)>0){
    print(paste('number of low outliers:', sum(data<=extreme.threshold.lower, na.rm = T)))
    }
    if (sum(data>=extreme.threshold.upper, na.rm = T)){
    print(paste('number of high outliers:', sum(data>=extreme.threshold.upper, na.rm = T)))
    }
    data[which(data<=extreme.threshold.lower)] = NA
    data[which(data>=extreme.threshold.upper)] = NA
    selDataSetOutliersRem[,varName] = data
  }
  return(selDataSetOutliersRem)
}