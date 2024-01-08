#!/bin/bash

inodes="i[001-022]"
ncpus=8
mem="8G"
time="2:00:00"
queue="tinyq"
njobs=$ncpus

rdir="results/ML_final"

export PYTHONHASHSEED=0
export SLURM_ARRAY_TASK_ID=1

for models in "ML_final_V1_corrected_combat_groupByExactDate_thr_V1_V1_LR_L2" "ML_final_V1_corrected_combat_groupByExactDate_thr_permuteLabels_V1_V1_LR_L2" "ML_final_V1_corrected_combat_groupByExactDate_thr_subsample0.9_V1_V1_LR_L2" "ML_final_V1_corrected_combat_groupByExactDate_thr_subsample0.9_permuteLabels_V1_V1_LR_L2"
do
      echo $models
      JOB_ID=1
      n=`grep . $rdir/$models -c`
      echo "There are $n models $rdir/$models"

      JOB_ID=$(sbatch --parsable --array=1-${n}%1 --wrap "python ML.py $models $njobs" \
      --job-name=${models} --output=logs/ML/${models}_%a.log \
      --mem=$mem --cpus-per-task=$ncpus --time=$time --partition=$queue --qos=$queue)

      sbatch --dependency=afterok:${JOB_ID} --kill-on-invalid-dep=yes --wrap "python ML_collect.py $models" \
     --job-name=collect_${models} --output=logs/collect_${models}.log \
     --mem=$mem --cpus-per-task=$ncpus --time=$time --partition=$queue --qos=$queue
done

LINE_N=13
models="ML_final_V1_corrected_combat_groupByExactDate_thr_V1_V1_LR_L2"
head -n ${LINE_N} $rdir/$models | tail -n 1
sbatch --array=${LINE_N} --wrap "python ML.py $models $njobs --save_data" \
--job-name="save_ML_data" --output=logs/save_ML_data.log \
--mem=$mem --cpus-per-task=$ncpus --time=$time --partition=$queue --qos=$queue
