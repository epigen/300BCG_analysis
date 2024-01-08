# 300BCG metadata

Metadata needed to run some of the scripts in this GitHub repository.

`complete_metadata.corrected.csv`

This is the main metadata file, which contains all sample annotations as well as personal immune profiles processed in different ways required for different types of analyses (for a clean copy of the personal immune profiles, please visit [Zenodo](https://doi.org/10.5281/zenodo.10288920)). A brief explanation regarding the metadata column names:

  - **CM** = circulating inflammatory markers
  - **CYTO** = cytokines (the quality of the cytokine data is suffixed as "\_good", "\_ok", "\_poor", "\_bad" and "_excluded"; we used only variables with the "\_good" quality suffix in the manuscript, but all are provided for completeness)
  - **PBMC\_\*** = cell type percentages or cell counts in PBMCs
  - **WB\_\*** = cell type percentages or cell counts in whole blood
  - **\*\_PERC** = cell type percentages
  - **\*\_PER\_ML** = cell counts (million cells per ml)
  - **\*\_SYSMEX\_PERC** = cell type percentages measured by Sysmex hematology analyzer (XN-450)
  - **\*\_SYSMEX\_PER\_ML** = cell counts (million cells per ml) measured by Sysmex hematology analyzer (XN-450)
  - **thm.\*** = donor scores ("thm" refers to "top half mean")
  - **thm.\*_V3_FC1.2_responder** = binarized (responder vs. non-responder) donor scores using a day-90 (V3) fold-change cut-off of 1.2
  - **IC\_\*** = value at the time of inoculation (day 0, V1)
  - **LFC\_\***, **LFC\_V2\_\***, **LFC\_V3\_\*** = log2 fold-changes, and log2 fold-changes at day 14 (V2) and day 90 (V3), respectively
  - **CORR\_\*** = corrected (e.g. corrected for sex, age, and blood-draw time)
  - **V1\_CORR\_\*** = corrected only using day 0 values (e.g. corrected for blood-draw time at day 0)

`gene_set_libraries`

This folder contains the gene/region sets used for enrichment analyses.
