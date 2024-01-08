#!/bin/bash

peak_annot="../data/DE/peaks_filtered_PBMC.csv.gz"
#peak_annot="metadata/gene_set_libraries/peaks_filtered_PBMC.GENE_AND_DISTAL_10kb.csv.gz"
#peak_annot="metadata/gene_set_libraries/peaks_filtered_PBMC.TSS_PROXIMAL.csv.gz"
min_size=15
max_size=500
for lib in 'KEGG_2019_Human' 'GO_Biological_Process_2018'
do
  echo ${lib}
  python gmt_genes_to_regions.py \
  --gmt_in metadata/gene_set_libraries/${lib}.gmt \
  --gmt_out metadata/gene_set_libraries/${lib}_min${min_size}_max${max_size}_regions.gmt \
  --peak_annot ${peak_annot} \
  --min_size ${min_size} \
  --max_size ${max_size} \
  --save_filtered_gmt metadata/gene_set_libraries/${lib}_min${min_size}_max${max_size}.gmt
done

