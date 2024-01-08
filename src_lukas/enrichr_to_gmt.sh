#!/bin/bash

for lib in 'GO_Biological_Process_2018' 'KEGG_2019_Human'
do
  echo $lib
  python enrichr_to_gmt.py --fn metadata/gene_set_libraries/$lib --out metadata/gene_set_libraries/$lib.gmt --make_safe
done
