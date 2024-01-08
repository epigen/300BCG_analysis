#!/bin/bash

cd ~/work/resources/LOLA/LOLARoadmap/hg19/Roadmap_Epigenomics_r9/regions
for x in *.narrowPeak
do
  ~/work/tools/liftOver/liftOver $x ~/work/tools/liftOver/hg19ToHg38.over.chain.gz ../../../hg38/Roadmap_Epigenomics_r9/regions/$x ../../../hg38/Roadmap_Epigenomics_r9/regions/$x.unmapped
done
wc -l ../../../hg38/Roadmap_Epigenomics_r9/regions/*.unmapped
rm ../../../hg38/Roadmap_Epigenomics_r9/regions/*.unmapped
