# 300BCG metadata

Metadata needed to run some of the scripts in this GitHub repository.


##complete_metadata.corrected.csv
This is the main file, containing most metadata. Many of the column names are self-explanatory, but a few important notes:

- The following abbraviations are used in the column names:
 - CM = circulating mediators
 - CYTO = cytokines
 - WB\_PER\_ML = Whole blood cell counts
 - SYSMEX\_PERC = Percentages of cells as measured by Sysmex
 - PBMC\_PERC = Percentages of cells in the PBMC fraction
- The quality of the cytokine data is suffixed as "\_good", "\_ok", "\_poor", "\_bad" and "_excluded". We focussed on the ones with "\_good" quality in the manuscript, but all are provided for completeness.
- The prefix "IC"
- The prefix "LFC" is an abbreciation of log fold-change, and "V2"/"V3" indicates the visit number (V2=day 14 and V3=day 90)

##gene\_set\_libraries
This folder contains the gene/region sets used in the paper.