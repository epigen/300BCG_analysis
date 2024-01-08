#!/bin/bash

ncpus=16
mem="256G"
time="12:00:00"
partition="shortq"

data_fn='../data/scRNAseq_YangLi/bcg0712_complete.QC2.SCT_log_counts.UMAP.QC3.celltypes.h5ad'
results_dir='results/scRNAseq/with_log_counts'

mkdir -p ${results_dir}

celltype_col='PBMC'
for reference in 'T0_RPMI' 'T0_LPS'
do
   job="scLM.${reference}.${celltype_col}"
   sbatch --wrap "python DE_scRNAseq.py --data_fn ${data_fn} --celltype_col ${celltype_col} --reference ${reference} --results_dir ${results_dir}" \
   --job-name=${job} --output=${results_dir}/${job}.log \
   --mem=${mem} --cpus-per-task=${ncpus} --time=${time} --partition=${partition} --qos=${partition}
done

for reference in 'interaction'
do
   job="scLM.${reference}.${celltype_col}"
   sbatch --wrap "python DE_scRNAseq.py --data_fn ${data_fn} --celltype_col ${celltype_col} --interaction_model --results_dir ${results_dir}" \
   --job-name=${job} --output=${results_dir}/${job}.log \
   --mem=${mem} --cpus-per-task=${ncpus} --time=${time} --partition=${partition} --qos=${partition}
done

celltype_col='celltypist'
for reference in 'T0_RPMI' 'T0_LPS'
do
    job="scLM.${reference}.${celltype_col}.subsample"
    sbatch --wrap "python DE_scRNAseq.py --data_fn ${data_fn} --celltype_col ${celltype_col} --reference ${reference} --subsample_celltypes --results_dir ${results_dir}" \
    --job-name=${job} --output=${results_dir}/${job}.log \
    --mem=${mem} --cpus-per-task=${ncpus} --time=${time} --partition=${partition} --qos=${partition}
done
