#!/bin/bash

vcf_dir="$1"
chr="$2"

picard_dir="../tools/picard_2.23.4"
picard="${picard_dir}/picard.jar"
chain="${picard_dir}/hg19ToHg38.over.chain.gz"
reference="${picard_dir}/hg38/hg38.fa"

input="${vcf_dir}/chr_${chr}_annotated.vcf.gz"
chr_input="${vcf_dir}/chr_${chr}_annotated_with_chr.vcf"
output="${vcf_dir}/GRCh38_chr_${chr}_annotated.vcf.gz"
reject="${vcf_dir}/GRCh38_chr_${chr}_annotated.vcf.rejected.gz"

zcat $input | awk '{if($0 !~ /^#/) print "chr"$0; else print $0}' - > ${chr_input}
java -Xmx120G -jar $picard LiftoverVcf I=${chr_input} O=$output CHAIN=$chain REJECT=$reject R=$reference WARN_ON_MISSING_CONTIG=true
rm ${chr_input}
