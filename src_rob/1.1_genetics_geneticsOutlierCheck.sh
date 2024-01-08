#!/bin/bash

#NOTE1: Unfortunately we were not allowed to share the genetics data due to privacy concerns. This code is provided for reference, but cannot be reproducibly run on the exact data used in the study. You can use your own or other public data.
#NOTE2: many of these steps will run most efficient if parallelized on a compute cluster 

#Set some paths
x1000GenomesFolder=""# set folder to get 1000genomes data, e.g. /nobackup/lab_bock/resources/1000G
x300bcgGeneticsFolder=""# set folder with 300BCG genetics data, e.g.  /research/lab_bock/projects/300BCG/data/SNP

##########PREPARE 1000G files##########

cd ${x1000GenomesFolder}

#1, Download the files as VCF.gz (and tab-indices)
prefix="ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chr" ;
suffix=".phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz" ;

for chr in {1..22}; do
    wget "${prefix}""${chr}""${suffix}" "${prefix}""${chr}""${suffix}".tbi ;
done

#2, Download 1000 Genomes PED file
wget ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/working/20130606_sample_info/20130606_g1k.ped ;

#3, Download the GRCh37 / hg19 reference genome
wget http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/human_g1k_v37.fasta.gz ;
wget http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/human_g1k_v37.fasta.fai ;

gunzip human_g1k_v37.fasta.gz ;

#4, Convert the 1000 Genomes files to BCF
for chr in {1..22}; do
	bcftools norm -m-any --check-ref w -f human_g1k_v37.fasta \
    ALL.chr"${chr}".phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz | \
    bcftools annotate -x ID -I +'%CHROM:%POS:%REF:%ALT' | \
    bcftools norm -Ob --rm-dup both \
	> ALL.chr"${chr}".phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.bcf ;
	
	bcftools index ALL.chr"${chr}".phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.bcf ;
done

#5, Convert the BCF files to PLINK format
for chr in {1..22}; do
	plink --noweb \
      --bcf ALL.chr"${chr}".phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.bcf \
      --keep-allele-order \
      --vcf-idspace-to _ \
      --const-fid \
      --allow-extra-chr 0 \
      --split-x b37 no-fail \
      --make-bed \
      --out ALL.chr"${chr}".phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes ;
done

#6, Prune variants from each chromosome
mkdir Pruned ;
for chr in {1..22}; do
	plink --noweb \
      --bfile ALL.chr"${chr}".phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes \
      --maf 0.10 --indep 50 5 1.5 \
      --out Pruned/ALL.chr"${chr}".phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes ;

	plink --noweb \
      --bfile ALL.chr"${chr}".phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes \
      --extract Pruned/ALL.chr"${chr}".phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.prune.in \
      --make-bed \
      --out Pruned/ALL.chr"${chr}".phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes ;
done

#7, Get a list of all PLINK files
find . -name "*.bim" | grep -e "Pruned" > ForMerge.list ;
sed -i 's/.bim//g' ForMerge.list ;

#8, Merge all projects into a single PLINK file
plink --merge-list ForMerge.list --out Merge ;




##########PREPARE 300BCG files##########
#This is done in a very similar way
#1, Convert to BCF

cd x300bcgGeneticsFolder
for chr in {1..22}; do
	bcftools norm -m-any --check-ref w -f ${x1000GenomesFolder}/human_g1k_v37.fasta \
    300BCG_chr_${chr}_annotated.vcf.gz | \
    bcftools annotate -x ID -I +'%CHROM:%POS:%REF:%ALT' | \
    bcftools norm -Ob --rm-dup both \
	> ${x300bcgGeneticsFolder}/300BCG_chr_${chr}_annotated.bcf ;

	bcftools index ${x300bcgGeneticsFolder}/300BCG_chr_${chr}_annotated.bcf 
done

#2, Convert to plink
for chr in {1..22}; do
	plink --noweb \
      --bcf 300BCG_chr_${chr}_annotated.bcf \
      --keep-allele-order \
      --vcf-idspace-to _ \
      --const-fid \
      --allow-extra-chr 0 \
      --split-x b37 no-fail \
      --make-bed \
      --out 300BCG_chr_${chr}_annotated ;
done

#3, Prune
for chr in {1..22}; do
    plink --noweb \
      --bfile 300BCG_chr_${chr}_annotated \
      --maf 0.10 --indep 50 5 1.5 \
      --out Pruned/300BCG_chr_${chr}_annotated ;


    plink --noweb \
      --bfile 300BCG_chr_${chr}_annotated \
      --extract Pruned/300BCG_chr_${chr}_annotated.prune.in \
      --make-bed \
      --out Pruned/300BCG_chr_${chr}_annotated.pruned ;
done

#4, Merge
find . -name "*.bim" | grep -e "Pruned" > ForMerge.list ;
sed -i 's/.bim//g' ForMerge.list ;
plink --merge-list ForMerge.list --out Merge ;




#######Merge both datasets########

#Some prep work
#1, Overlap common SNPS
##——RUN THIS SIMPLE CODE IN R!--
#library(data.table)
#map2 = fread("/scratch/lab_bock/rterhorst/Research/Resources/300bcg_genetics/Merge.bim", header=F)
#map1 = fread("/scratch/lab_bock/rterhorst/Research/Resources/1000G/Merge.bim", header=F)
#common.snps = which(map2$V2 %in% map1$V2)
#write.table(map2$V2[common.snps], file="/scratch/lab_bock/rterhorst/Research/Resources/300bcg_genetics/commonSnps_1000G_300Bcg.snps", sep="\t", col.names=F, row.names=F, quote=F )
##——END OF CODE IN R--

#2, Merge with 300BCG data
file1000G=${x1000GenomesFolder}"/Merge"
file300Bcg=${x300bcgGeneticsFolder}"/Merge"
listSnps=${x300bcgGeneticsFolder}"/commonSnps_1000G_300Bcg.snps"
fileOut=${x300bcgGeneticsFolder}"/merged_1000G_300Bcg"

plink --bfile ${file1000G} --bmerge ${file300Bcg}.bed ${file300Bcg}.bim ${file300Bcg}.fam --extract ${listSnps} --make-bed -out ${fileOut}

#3. Check relatedness
plink --bfile ${file300Bcg} --genome



##########Make plots##########
#If you want to run PCA and MDS yourself
plink --bfile ${fileOut} --pca
plink --bfile ${fileOut} --cluster --mds-plot 10

#Or use the R script below
#1.2_geneticsEthnicityPlots.R contains more detailed plot info
