# 300BCG code src\_lukas

Description of scripts to run in order to create some of the figures from the manuscript
entitled **"Multi-omics analysis of innate and adaptive responses to BCG vaccination reveals epigenetic cell states that predict trained immunity"** by Moorlag, Folkman, ter Horst, Krausgruber et al., *Immunity* (2024). doi: https://doi.org/10.1016/j.immuni.2023.12.005. This only concerns the code stored in `src_lukas`.

Note: The analyses for this paper were performed by Rob ter Horst (`src_rob`) and Lukas Folkman (`src_lukas`).

Specifically, the following figures can be reproduced (before formatting) using the code in `src_lukas`:

- Figures 1B–E
- Figures 2B–J
- Figures 4B–I
- Figures 5B–D
- Figures 6B–E
- Figures S1A–D
- Figures S2A–F
- Figures S4A–D
- Figures S5A–B and S5D–G

## Part 1: Statistical analysis of immune cells, inflammatory markers, and cytokines

NOTE: All `job_*.sh` scripts are written for SLURM workload manager

#### Dependencies

    conda env create -f envs/300BCG_lukas.yml

#### Run analyses (SLURM jobs) with linear models and linear mixed models

    conda activate 300BCG_lukas
    bash job_LM_var.sh
    bash job_LM_day0_factors.sh
    bash job_LM_day90_factors.sh
    bash job_LM_seasonal.sh
    bash job_LM_trained.sh

#### Summarize results, make plots and tables

    jupyter notebook Figures_immune_overview.ipynb
    jupyter notebook Figures_immune_day0.ipynb
    jupyter notebook Figures_immune_trained.ipynb

## Part 2: Statistical analysis of the chromatin accessibility profiles

NOTE: All `job_*.sh` scripts are written for SLURM workload manager

#### Dependencies

    conda env create -f envs/300BCG_lukas.yml

#### Select features for downstream analysis

    conda activate 300BCG_lukas
    python data_select_features.py

#### Prepare gene/region set databases (KEGG and GO Biological Process, TFBS from HOCOMOCO)

    conda activate 300BCG_lukas
    bash enrichr_to_gmt.sh
    bash gmt_genes_to_regions.sh
    python select_regions_mapped_to_genes.py
    bash job_TFBS_FIMO.sh

#### Remap genomic regions and genes from publicly available studies for validation

    conda activate 300BCG_lukas
    python remap_previous_studies.py

#### Make configuration files for all differential accessibility analyses

    jupyter notebook LIMMA_config.ipynb

#### Run analyses with LIMMA, perform gene/region set enrichments

    conda activate 300BCG_lukas
    bash job_DREAM_var.sh
    bash job_LIMMA_day0_factors.sh
    bash job_LIMMA_day0_immune.sh
    bash job_LIMMA_day90_immune.sh
    bash job_LIMMA_seasonal.sh
    bash job_LIMMA_trained.sh
    bash job_LIMMA_trained_bootstrap.sh

#### Summarize results, make plots and tables

    jupyter notebook Figures_ATACseq_overview.ipynb
    jupyter notebook Figures_ATACseq_day0.ipynb
    jupyter notebook Figures_ATACseq_day0_suppl_table.ipynb
    jupyter notebook Figures_ATACseq_seasonal.ipynb
    jupyter notebook Figures_ATACseq_trained.ipynb

#### Prepare, run, and summarize the machine-learning analysis

    conda activate 300BCG_lukas
    python ML_data.py
    python ML_config.py
    bash job_ML.sh
    jupyter notebook Figures_ML.ipynb

## Part 3: Statistical analysis of the scRNA-seq profiles

NOTE: All `job_*.sh` scripts are written for SLURM workload manager

#### Dependencies

    conda env create -f envs/300BCG_lukas.yml
    conda env create -f envs/sctransform.yml
    conda env create -f envs/celltypist.yml

#### Preprocess scRNA-seq data

    jupyter notebook scRNAseq_1_read_RDS.ipynb
    jupyter notebook scRNAseq_2_preprocess.ipynb
    jupyter notebook scRNAseq_3_SCT.ipynb
    jupyter notebook scRNAseq_4_read_SCT.ipynb
    jupyter notebook scRNAseq_5_unsupervised.ipynb
    jupyter notebook scRNAseq_6_celltypes.ipynb

#### Run differential expression analysis

    conda activate 300BCG_lukas
    bash job_DE_scRNAseq.sh

#### Summarize results, make plots and tables

    jupyter notebook Figures_scRNAseq.ipynb
    jupyter notebook Figures_scRNAseq_suppl.ipynb
    jupyter notebook Figures_scRNAseq_suppl_table.ipynb
