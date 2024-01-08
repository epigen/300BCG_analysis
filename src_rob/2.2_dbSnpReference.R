library("VariantAnnotation")
library("vcfR")
library(snpStats)
library(R.utils)
library(data.table)
library(stringr)
library(dplyr)

source('processingFxnsGeneral.R')
source('plottingFxnsGeneral.R')
source('300bcg_generalFxns.R')


###Get the variables we need
library(config)
Sys.setenv(R_CONFIG_ACTIVE = "300bcg")
config <- config::get()
for (x in names(config)){eval(parse(text=paste0(x,"='",config[[x]],"'")))}

###(1) Run the following commands on the Unix command line! Requires bcftools
#cd /Volumes/nobackup/lab_bock/projects/300BCG/resources/dbSNP/ ##REPLACE WITH YOUR FOLDER, i.e. resourceDir/dbSNP/
#wget https://ftp.ncbi.nih.gov/snp/organisms/human_9606_b151_GRCh38p7/VCF/00-common_all.vcf.gz
#wget https://ftp.ncbi.nih.gov/snp/organisms/human_9606_b151_GRCh38p7/VCF/00-common_all.vcf.gz.tbi
#bcftools view --no-header 00-common_all.vcf.gz | cut -f 1-5 > human_9606_b151_GRCh38p7_SNPdb.tdf

###(2) Process files
###Load new RS-id mapping
print('loading dbSNP')
rsInfo = fread(file.path(resourceDir,"dbSNP","human_9606_b151_GRCh38p7_SNPdb.tdf"))
##Keep just the SNPs and remove indels
print('selecting rows')
rowsToKeep = which(unlist(lapply(rsInfo[,4], nchar))==1 & unlist(lapply(rsInfo[,5], nchar))==1)
rsInfo = rsInfo[rowsToKeep,]
colnames(rsInfo) = c('CHROM','POS','ID','REF','ALT')
head(rsInfo,20)
print('writing file')
fwrite(x = rsInfo, file = file.path(resourceDir,"dbSNP","human_9606_b151_GRCh38p7_SNPdb.filt.tdf"), col.names = T, row.names = F)

for (chrSel in c(1:22)){
  print (paste0('chr=',chrSel))
  rsInfo.subset = rsInfo[J(as.character(chrSel)), nomatch=0L, on = .(CHROM)]
  fwrite(x = rsInfo.subset, file = file.path(resourceDir,"dbSNP",paste0("human_9606_b151_GRCh38p7_SNPdb.filt_chr",chrSel,".tdf")), col.names = T, row.names = F)
}