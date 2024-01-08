# 300BCG code src\_rob

Scripts to run in order to create the basic figures from the manuscript
entitled **"Multi-omics analysis of innate and adaptive responses to BCG vaccination reveals epigenetic cell states that predict trained immunity"** by Moorlag, Folkman, ter Horst, Krausgruber et al., *Immunity* (2024). doi: https://doi.org/10.1016/j.immuni.2023.12.005. This only concerns the code stored in `src_rob`.

Note: The analyses for this paper were performed by Rob ter Horst (`src_rob`) and Lukas Folkman (`src_lukas`).

Specifically, the following figures can be reproduced (before formatting) using the code in `src_rob`:

- Figure 1F
- Figure 1G
- Figure 1H (just plotting, analysis in `src_lukas`)
- Figure 1I (just plotting, analysis in `src_lukas`)
- Figure 1J
- Figure 3B
- Figure 3C
- Figure 3D
- Figure 3E (just plotting, analysis in `src_lukas`)
- Figure 3F (just plotting, analysis in `src_lukas`)
- Figure 3G (just plotting, analysis in `src_lukas`)
- Figure 3H
- Figure 5E
- Figure 5F
- Figure S1E
- Figure S1F
- Figure S1G (just plotting, analysis in `src_lukas`)
- Figure S1H (just plotting, analysis in `src_lukas`)
- Figure S1I
- Figure S1J
- Figure S3A
- Figure S3B (just plotting, analysis in `src_lukas`)
- Figure S3C
- Figure S5C

The .yml file `envs/300BCG_rob.yml` provides an overview of the packages installed in the conda environment used to run this code. You can use this file to create a similar conda environment on your machine.

NOTE: All `*.batch` scripts are written for SLURM workload manager, which in some cases were designed to be submitted by an ".R" or a ".submit" script. Below the scripts are listed that can be used for each step in the analysis. If you do not use SLURM, some of the scripts might have to be slightly adapted to run on your system.

NOTE2: Please change the paths in the file "config.yml"

# **GENETIC DATA**
## 1. Outlier testing genetics
NOTE: Outlier testing is done on the GrCH37 data, since the 1000G is also in this version

- Step 1. The data processing is done with: `1.1_genetics_geneticsOutlierCheck.sh`
- Step 2. The ethnicity check plots can be done with: `1.2_geneticsEthnicityPlots.R`

## 2. Processing the genetics
- Step1. Lifting the genetics from GrCH37 to GrCH38 --> see src\_lukas
- Step2. Filter the data for MAF HWE and the outliers: `2.1_filterGeneticsData_vcf.sh`
- Step3. Download and merge the latest mapping of rs-id to postion for GrCH38: `2.2_dbSnpReference.R`
- Step4. Convert the genetics data to dosage format and save: `2.3_convertVcfToDosages.R` / `2.3_convertVcfToDosages.batch` / `2.3_convertVcfToDosages.submit`

## 3. cQTL analyses
- Step1. Run the QTL analysis: `3.1_QTL.R` / `3.1_QTL.batch`
- Step2. Maps the QTL sets to .gmt files with entrez ids instead of gene names, needed for some of the enrichment analyses: `3.2_QTLEnrichment_mapGenesymbolToEntrez.R`
- Step3. `3.3_QTLEnrichment.R` / `3.3_QTLEnrichment.batch` - QTL enrichment PASCAL --> runs both the mapping to top genes and the actual GSEA [before running this I needed to convert our custom sets to entrez using `3.2_QTLEnrichment_mapGenesymbolToEntrez.R`]

## 4. Plots genetics
- manhattan plots & top hit example plots [Tables & Figures]: `4.1_topCQtlHits.R` / `4.1_topCQtlHits.batch`
- pathway bar plots [Tables & Figures]: `4.2_QTLEnrichment_combineResults.R` / `4.2_QTLEnrichment_combineResults.batch`


# **EPIGENETIC DATA**
## 5. Plots epigenetics
- box-plots: `5.1_plotTopHitsAtac.R` / `5.1_plotTopHitsAtac.batch`
- example plots showing the top hits: `5.2_examplePlotsAtacCytoRegres.R`
- plot the pathway bar plots: `5.3_plotTopPathwaysEpigeneticsLukas.R` / `5.3_plotTopPathwaysEpigeneticsLukas.batch`

# **GENETIC + EPIGENETIC DATA**
## 6. Bayesian
- Step 1. Bayesian percentage of variance calculations: `6.1_atacGeneticsBayesian.R` / `6.1_atacGeneticsBayesian.batch` / `6.1_atacGeneticsBayesian.submit`
- Step 2. Bayesian percentage of variance plots: `6.2_atacGeneticsBayesian_plots.R` /`6.2_atacGeneticsBayesian_plots.batch` [Tables & Figures]

## 7. Overlap genetics and epigenetics
- Step1. create the haplotype blocks: `7.1_createHaplotypeBlocks.R` / `7.1_createHaplotypeBlocks.submit` / `7.1_createHaplotypeBlocks.batch`
- Step2. extend regions with the haplotype blocks: `7.2_expandPeaksWithHaplotypeBlocks.R`
- Step3. calculate the overlap: `7.3_overlapGeneticsEpigenetics.R` / `7.3_overlapGeneticsEpigenetics.batch`  [Tables & Figures]

# **VALIDATION USING PUBLIC TRAINED IMMUNITY DATA**
## 8. Public trained immunity overlap
### Calculations
- Step 1. Overlap our ATAC-seq results with published gene/regions sets and scores linked to trained immunity: `8.1_overlapAtacWithPublicData.R` [check associated website for the region sets]

### Plots
- Step 2. Heatmap plots: `8.2_overlapAtacWithPublicData_plots.R`
- Step 3. Example GSEA plots: `8.3_overlapAtacWithPublicData_plotGseaExamples.R`

# GENERAL FUNCTIONS
Some general functions that can be sourced from the following scripts:

 - `plottingFxnsGeneral.R`
 - `300bcg_generalFxns.R`
 - `processingFxnsGeneral.R`
