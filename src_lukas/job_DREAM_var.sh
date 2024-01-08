#!/bin/bash

snodes="s[001-004]"
bnodes="b[001-003],b[006-007],b[011-012],b[014-018]"
ncpus=8
mem="8G"
time="2:00:00"
partition="tinyq"

celltype="PBMC"
results="results_check"

for model in 'varPartTop.batch_TSS_enr_blood_corrected.donor'
do
  JOB_ID=1
  ls ${results}/DE/${celltype}.${model}/de_config_${celltype}.${model}.yml

# DREAM
  job="varPart_${celltype}.${model}"
  printf "${job}: "
  JOB_ID=$(sbatch --parsable --wrap "Rscript limma.R ${results}/DE/${celltype}.${model}/de_config_${celltype}.${model}.yml $ncpus" \
  --job-name=${job} --output=logs/${job}.log \
  --mem=16G --cpus-per-task=$ncpus --time=2-00:00:00 --partition=mediumq --qos=mediumq --exclude=$bnodes)
  printf "Submitted batch job ${JOB_ID}\n"

# ENRICHR
  job="Enrichr_${model}"
  sbatch --dependency=afterok:${JOB_ID} --kill-on-invalid-dep=yes \
  --wrap "python Enrichr_donorVar.py ${results}/DE/${celltype}.${model}/de_config_${celltype}.${model}.yml" \
  --job-name=${job} --output=logs/${job}.log \
  --mem=$mem --cpus-per-task=$ncpus --time=$time --partition=$partition --qos=$partition --exclude=$bnodes

# LOLA
  job="LOLA_${model}"
  sbatch --dependency=afterok:${JOB_ID} --kill-on-invalid-dep=yes \
  --wrap "python LOLA_donorVar.py ${results}/DE/${celltype}.${model}/de_config_${celltype}.${model}.yml" \
  --job-name=${job} --output=logs/${job}.log \
  --mem=180G --cpus-per-task=$ncpus --time=$time --partition=$partition --qos=$partition --exclude=$bnodes
done
