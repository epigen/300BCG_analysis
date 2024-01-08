#!/bin/bash

##Make sure you have tabix, htslib, bcftools and vcftools installed.
##You can also create a conda environment that contains these
conda activate 300bcg_main #activate env
#Install the packages below if not part of your cluster
module load tabix
module load htslib/1.2.1
module load bcftools
#module load vcftools #Now loaded from conda 

##Define where folder the vcf files are located
vcfFileDir=''#'/home/lfolkman/bcg/data/genotype/vcf'
cd ${vcfFileDir}

# I will use this code to filter the GRCh38 data for MAF and HWE AND remove the genetic outliers
folderToSave=mafHweOutlierFilteredGrCh38
mkdir ${folderToSave}

individualsToRemove=(6_6 40_40 41_41 58_58 74_74 106_106 139_139 196_196 200_200 201_201 208_208)
printf "%s\n" "${individualsToRemove[@]}" > indivToRemove.txt

for chrSel in {1..22}
do
	inputFile=GRCh38_chr_${chrSel}_annotated.vcf.gz
	outputFile=${folderToSave}/GRCh38_300BCG_chr_${chrSel}_annotated_MAF0.1_HWE1e-5_noOutliers
	## Do the filtering
	vcftools --gzvcf ${inputFile} --remove indivToRemove.txt --maf 0.1 --hwe 1e-5 --recode --out ${outputFile}
	##Update tags
	bcftools +fill-tags ${outputFile}.recode.vcf > ${outputFile}.recode.tags.vcf

	bgzip -c ${outputFile}.recode.vcf > ${outputFile}.recode.tags.vcf.gz
	tabix -p vcf ${outputFile}.recode.tags.vcf.gz
done

#####Some checks, not necessary to run
##Check which SNPs have actually been removed
#bcftools isec -p ${folderToSave}/isec_results_chr${chrSel} ${inputFile} ${outputFile}.recode.tags.vcf.gz
#
##Sanity check to see if the ones that were removed indeed had a MAF<.1 (or HWE problem)
#bcftools +fill-tags ${folderToSave}/isec_results_chr22/0000.vcf  -- -t MAF | tail -5
#bcftools +fill-tags ${folderToSave}/isec_results_chr22/0000.vcf | tail -5
