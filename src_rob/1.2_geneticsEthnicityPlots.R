library(openxlsx)
library(data.table)
library("R.utils")
library(ggplot2)
library(stringr)
print('starting')

###Get the variables we need
library(config)
Sys.setenv(R_CONFIG_ACTIVE = "300bcg")
config <- config::get()
for (x in names(config)){eval(parse(text=paste0(x,"='",config[[x]],"'")))}

###Define file locations of the files produced and downloaded with "1.1_geneticsOutlierCheck.sh" in config.yml

###Define which data to plot
#"merged_1000G_300Bcg" = 300bcg + 1000G
#"merge" = just 300 bcg

dataSel = "merged_1000G_300Bcg" 
selectedDimRedMethod = 'PCA' #either "MDS" or "PCA"
dim1Sel = 'PC1'
dim2Sel = 'PC2'

# read in the PED data
PED <- fread(file.path(x1000GenomesFolder,'20130606_g1k.ped'), header = TRUE, data.table = F)

#Plink data
library(snpStats)
geneticsData = read.plink(file.path(x300bcgGeneticsFolder,dataSel))

#mapping between pop and superpop
populationInfo = openxlsx::read.xlsx(xlsxFile = file.path('populationSuperPopulation.xlsx'))
mapPopToSuperPop = c(populationInfo$Super.Population.Code,'300BCG')
names(mapPopToSuperPop) = c(populationInfo$Population.Code,'300BCG')
mapPopToSuperPop[['FIN']] = 'FIN'
PED$superPopulation = unlist(lapply(PED$Population, function(x) mapPopToSuperPop[[x]]))

#Subset data if needed
oldwd = getwd()
setwd(x300bcgGeneticsFolder)

#Create a list of individuals to remove/keep
toRemove = c() #if we want to plot all

#Some specific groups we can remove:
#toRemoveEthnic = c("40_40", "41_41", "58_58", "106_106", "208_208")
#toRemove = c("6_6", "40_40", "41_41", "58_58", "74_74", "106_106", 
#             "139_139", "196_196", "200_200", "201_201", "208_208")

allIndiv = geneticsData$fam$member
toKeep = setdiff(allIndiv,toRemove)
write.table(x = cbind.data.frame('0',toKeep), file = 'toExtract.txt', sep = ' ', row.names = F, col.names = F, quote = F)  
system(paste0("/Users/rterhorst/local/plink_mac_20200616/plink --bfile ",dataSel," --keep toExtract.txt --make-bed --out ",dataSel,"_selected"))
selFile = paste0(dataSel,'_selected')

#perform PCA
system(paste0(plinkExecutable,' --bfile ',selFile,' --pca'))
system(paste0(plinkExecutable,' --bfile ',selFile,' --cluster --mds-plot 10'))
setwd(oldwd)

#Load the PCA result
#PCA or MDS
if (selectedDimRedMethod == 'PCA'){
  eigenvec <- read.table(file.path(x300bcgGeneticsFolder,'plink.eigenvec'), header = FALSE, skip=0, sep = ' ')
}else{
  eigenvec <- fread(file.path(x300bcgGeneticsFolder,'plink.mds'), header = TRUE, skip=0, data.table = F)
  eigenvec$SOL <- NULL
}

rownames(eigenvec) <- eigenvec[,2]
eigenvec <- eigenvec[,2:ncol(eigenvec)]
colnames(eigenvec)[2:ncol(eigenvec)] <- paste('PC', c(1:20), sep = '')

PED <- PED[which(PED$`Individual ID` %in% rownames(eigenvec)), ]

if (dataSel=='merge'){
  eigenvecWithPops = eigenvec
  eigenvecWithPops$`Individual ID`=as.character(eigenvecWithPops[,1])
  eigenvecWithPops$superPopulation = '300BCG'
  
  eigenvecWithPops$overallOutlier = NA
  
  topNr=6
  
  eigenvecWithPops$pc1Outlier = NA
  selPc = order(eigenvecWithPops[,dim1Sel], decreasing = T)[c(c(1:topNr),c((nrow(eigenvecWithPops)-topNr):nrow(eigenvecWithPops)))]
  eigenvecWithPops$pc1Outlier[selPc] = eigenvecWithPops$`Individual ID`[selPc]
  eigenvecWithPops$overallOutlier[selPc] = eigenvecWithPops$`Individual ID`[selPc]
  
  eigenvecWithPops$pc2Outlier = NA
  selPc = order(eigenvecWithPops[,dim2Sel], decreasing = T)[c(c(1:topNr),c((nrow(eigenvecWithPops)-topNr):nrow(eigenvecWithPops)))]
  eigenvecWithPops$pc2Outlier[selPc] = eigenvecWithPops$`Individual ID`[selPc]
  eigenvecWithPops$overallOutlier[selPc] = eigenvecWithPops$`Individual ID`[selPc]
  
}else{
  eigenvecWithPops = base::merge(PED[,c('Individual ID','Population')], eigenvec,by.x='Individual ID',by.y='V2',all=T)
  eigenvecWithPops$`Individual ID` = gsub('^[0-9]{1,3}\\_','300BCG_',eigenvecWithPops$`Individual ID`)
  
  eigenvecWithPops[startsWith(eigenvecWithPops$`Individual ID`,prefix = '300BCG'),'Population']='300BCG'
  
  
  eigenvecWithPops$superPopulation = unlist(lapply(eigenvecWithPops$Population, function(x) mapPopToSuperPop[[x]]))
  eigenvecWithPops$superPopulation = factor(eigenvecWithPops$superPopulation, levels=c(unique(mapPopToSuperPop)))
  eigenvecWithPops = eigenvecWithPops[rev(c(1:nrow(eigenvecWithPops))),,drop=F]
  
  eigenvecWithPops$overallOutlier = NA
  
  eigenvecWithPops$pc1Outlier = NA
  selPc = eigenvecWithPops[,dim1Sel]>-.005 & eigenvecWithPops$Population=='300BCG'
  eigenvecWithPops$pc1Outlier[selPc] = eigenvecWithPops$`Individual ID`[selPc]
  eigenvecWithPops$overallOutlier[selPc] = eigenvecWithPops$`Individual ID`[selPc]
  
  eigenvecWithPops$pc2Outlier = NA
  selPc = eigenvecWithPops[,dim2Sel]<.018 & eigenvecWithPops$Population=='300BCG'
  eigenvecWithPops$pc2Outlier[selPc] = eigenvecWithPops$`Individual ID`[selPc]
  eigenvecWithPops$overallOutlier[selPc] = eigenvecWithPops$`Individual ID`[selPc]
}

library(ggrepel)
b <- ggplot(eigenvecWithPops, aes_string(x = dim1Sel, y = dim2Sel)) + 
  geom_point(aes(color = superPopulation), alpha=.7) +
  geom_label_repel(aes(label = overallOutlier, size = 3.5))
print(b)