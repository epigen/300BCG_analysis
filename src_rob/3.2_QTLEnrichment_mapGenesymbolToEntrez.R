library(R.utils)
library(data.table)
library('org.Hs.eg.db')
library(DBI)
#I got some of the code from a blog post that gives a nice summary: https://www.gungorbudak.com/blog/2018/08/07/convert-gene-symbols-to-entrez-ids-in-r/

###Get the variables we need
library(config)
Sys.setenv(R_CONFIG_ACTIVE = "300bcg")
config <- config::get()
for (x in names(config)){eval(parse(text=paste0(x,"='",config[[x]],"'")))}

geneSetFolder = file.path(metaDataFolder,"gene_set_libraries")

# Names of the original files
geneSetFiles =  c("KEGG_2019_Human",  "GO_Biological_Process_2018")

#SYNOMYMOUS IDS
# use sql to get alias table and gene_info table (contains the symbols)
# first open the database connection
dbCon <- org.Hs.eg_dbconn()
# write your SQL query
sqlQuery <- 'SELECT * FROM alias, gene_info WHERE alias._id == gene_info._id;'
# execute the query on the database
aliasSymbol <- dbGetQuery(dbCon, sqlQuery)

for (baseFileName in geneSetFiles){
  entrezIdList= list()
  #Load gene set file
  geneSetSel = readLines(file.path(geneSetFolder,paste0(baseFileName,'.gmt')))
  rowIndex=0
  #iterate over each gene set and do the mapping
  for (rowSel in geneSetSel){
    rowIndex = rowIndex + 1
    print(paste('row',rowIndex,'out of',length(geneSetSel)))
    splitNames = strsplit(rowSel, split = '\t')[[1]]
    geneNames = splitNames[3:length(splitNames)]
    
    #For the orf gene names, the ORF part needs to be in lower case orf
    orfMatches = grep('.*[0-9]ORF[0-9].*',x = geneNames)
    if (length(orfMatches)>0){
      geneNames[orfMatches] = gsub('ORF','orf', x = geneNames[orfMatches])
    }
    
    entrezIds = mapIds(org.Hs.eg.db, geneNames, 'ENTREZID', 'SYMBOL', multiVals="list")
    #Unlist to named character vector
    entrezIds = unlist(entrezIds)
    print(entrezIds[order(names(entrezIds))])
    #There are some problems with some genes, e.g. RNR2 (or RNR1) vs RNR2-MT
    entrezIdsFull = entrezIds
    if (sum(is.na(entrezIds))>0){
      unmappedIds = names(entrezIds)[which(is.na(entrezIds))]
      #Make sure I have all aliases for the gene symbols
      newIDs = aliasSymbol[which(aliasSymbol[,'alias_symbol'] %in% unmappedIds),c('alias_symbol','symbol'), drop=F]
      if (nrow(newIDs)>0){
        #Map these symbols
        newMapping = mapIds(org.Hs.eg.db, newIDs$symbol, 'ENTREZID', 'SYMBOL')
        entrezIdsFull = c(entrezIds, newMapping)
        entrezIdsFull = entrezIdsFull[!is.na(entrezIdsFull)]
      }
    }
    
    entrezIdsFull = entrezIdsFull[!is.na(entrezIdsFull)]
    
    #filter for chrom 1-22 (that gets rid of any of the double mitochondrial ones as well)
    x <- org.Hs.egCHR
    chromOfEntrez = as.list(x[entrezIdsFull])
    toKeepIds = names(chromOfEntrez[chromOfEntrez %in% as.character(c(1:22))])
    
    print(setdiff(entrezIdsFull,toKeepIds))
    pathwayName = splitNames[[1]]
    descript = 'NoDescript'
    entrezIdList[[rowIndex]] = paste(pathwayName,descript, paste(entrezIds, collapse = '\t'), sep='\t')
  }
  entrezIdList = as.character(entrezIdList)
  
  #Write the new gene sets
  writeLines(text = entrezIdList, con = file.path(geneSetFolder,paste0(baseFileName,'_entrez.gmt')))
}

