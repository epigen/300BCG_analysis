#!/bin/bash

snodes="s[001-004]"
bnodes="b[001-003],b[006-007],b[011-012],b[014-018]"
ncpus=8
mem="8G"
time="2:00:00"
partition="tinyq"

celltype="PBMC"
results="results_check"

for model in "fixed.donor_as_mixed.batch_corrected.sex.age.blood.TSS_enr.visit_time.season_interaction"
do
  JOB_ID=1
  ls ${results}/DE/${celltype}.${model}/de_config_${celltype}.${model}.yml
     
#  LIMMA
  job="limma_${celltype}.${model}"
  printf "${job}: "
  JOB_ID=$(sbatch --parsable --wrap "Rscript limma.R ${results}/DE/${celltype}.${model}/de_config_${celltype}.${model}.yml $ncpus" \
  --job-name=${job} --output=logs/${job}.log \
  --mem=24G --cpus-per-task=$ncpus --time=2-00:00:00 --partition=mediumq --qos=mediumq --exclude=$bnodes)
  printf "Submitted batch job ${JOB_ID}\n"

#  VOLCANO
  job="volcano_${celltype}.${model}"
  printf "${job}: "
  JOB_ID=$(sbatch --parsable --dependency=afterok:${JOB_ID} --kill-on-invalid-dep=yes --wrap "python DE_volcanos.py --config ${results}/DE/${celltype}.${model}/de_config_${celltype}.${model}.yml" \
  --job-name=${job} --output=logs/${job}.log \
  --mem=$mem --cpus-per-task=$ncpus --time=$time --partition=$partition --qos=$partition)

#  COMBINE F-TEST RESULTS
  job="combine_${celltype}.${model}"
  printf "${job}: "
  JOB_ID=$(sbatch --parsable --dependency=afterok:${JOB_ID} --kill-on-invalid-dep=yes \
  --wrap "python DE_combine_season_F_test_results.py ${results}/DE/${celltype}.${model}/de_config_${celltype}.${model}.yml" \
  --job-name=${job} --output=logs/${job}.log \
  --mem=$mem --cpus-per-task=$ncpus --time=$time --partition=$partition --qos=$partition)

#  CLUSTERING
  job="seasonalClusters_${celltype}.${model}"
  JOB_ID=$(sbatch --parsable --dependency=afterok:${JOB_ID} --kill-on-invalid-dep=yes \
  --wrap "python DE_seasonalClusters.py ${results}/DE/${celltype}.${model}/de_config_${celltype}.${model}.yml" \
  --job-name=$job --output=logs/$job.log \
  --mem=$mem --cpus-per-task=$ncpus --time=$time --partition=$partition --qos=$partition --exclude=$bnodes)

#  ENRICHR
   job="Enrichr_seasonalClusters_${celltype}.${model}"
   sbatch --dependency=afterok:${JOB_ID} --kill-on-invalid-dep=yes \
   --wrap "python Enrichr_seasonalClusters.py ${results}/DE/${celltype}.${model}/de_config_${celltype}.${model}.yml" \
   --job-name=${job} --output=logs/${job}.log \
   --mem=$mem --cpus-per-task=$ncpus --time=$time --partition=$partition --qos=$partition

#   LOLA
   # I run this per cluster because otherwise it runs out of memory with the large HOCOMOCO database
   for visit in "V2" "V3"
   do
     for clust in 1 2 3 4 5
     do
       job="LOLA_seasonalClusters.${visit}.clust${clust}.${celltype}.${model}"
       sbatch --dependency=afterok:${JOB_ID} --kill-on-invalid-dep=yes \
       --wrap "python LOLA_seasonalClusters.py ${results}/DE/${celltype}.${model}/de_config_${celltype}.${model}.yml ${visit} ${clust}" \
       --job-name=${job} --output=logs/${job}.log \
       --mem=180G --cpus-per-task=$ncpus --time=$time --partition=$partition --qos=$partition --exclude=$bnodes
     done
   done

   echo ""
done
