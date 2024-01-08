#!/bin/bash

ncpus=8
mem="8G"
time="2:00:00"
partition="tinyq"

celltype="PBMC"
results="results_bootstrap"

for phenotype in "thm.innate_nonspecific_24h_wo_LAC_IL10_IL1ra_V3_FC1.2"
do
  JOB_ID=1

  for seed in `seq 101 1000`
  do

  model="BS${seed}.donor_as_mixed.batch.sex.age.blood.TSS_enr.visit_time.${phenotype}_responder"

  # LIMMA
    job="limma_${celltype}.${model}"
    JOB_ID=$(sbatch --parsable --wrap "Rscript limma.R ${results}/DE/${celltype}.${model}/de_config_${celltype}.${model}.yml ${ncpus}" \
    --job-name=${job} --output=logs/${job}.log \
    --mem=16G --cpus-per-task=$ncpus --time=2-00:00:00 --partition=mediumq --qos=mediumq)
    printf "Submitted batch job ${JOB_ID}\n"

  # VOLCANO
    job="volcano_${celltype}.${model}"
    printf "${job}: "
    JOB_ID=$(sbatch --parsable --dependency=afterok:${JOB_ID} --kill-on-invalid-dep=yes --wrap "python DE_volcanos.py --config ${results}/DE/${celltype}.${model}/de_config_${celltype}.${model}.yml" \
    --job-name=${job} --output=logs/${job}.log \
    --mem=$mem --cpus-per-task=$ncpus --time=$time --partition=$partition --qos=$partition)

  # ENRICHR
    job="Enrichr_${celltype}.${model}"
    printf "${job}: "
    sbatch --dependency=afterok:${JOB_ID} --kill-on-invalid-dep=yes --wrap "python DE_enrichr.py --config ${results}/DE/${celltype}.${model}/de_config_${celltype}.${model}.yml" \
    --job-name=${job} --output=logs/${job}.log \
    --mem=$mem --cpus-per-task=$ncpus --time=$time --partition=$partition --qos=$partition

  #  # LOLA
  #  lola_universe=resources/LOLA/lola_universe.bed  # resources/LOLA/lola_universe.TSS_PROXIMAL.bed
  #  region_filter=""  # "--regions TSS_PROXIMAL"
  #  job="LOLA_${celltype}.${model}"
  #  echo ${job}
  #  echo "Universe size:"
  #  wc -l ${lola_universe}
  #  sbatch --dependency=afterok:${JOB_ID} --kill-on-invalid-dep=yes \
  #  --wrap "python DE_lola.py --top_n 1000 --config ${results}/DE/${celltype}.${model}/de_config_${celltype}.${model}.yml --n_jobs ${ncpus} ${region_filter} --universe ${lola_universe}" \
  #  --job-name=${job} --output=logs/${job}.log \
  #  --mem=180G --cpus-per-task=$ncpus --time=$time --partition=$partition --qos=$partition

  done
done
