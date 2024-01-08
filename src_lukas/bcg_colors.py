import numpy as np

BLUE = '#4c72b0'
ORANGE = '#dd8452'
GREEN = '#55a868'
RED = '#c44e52'
PURPLE = '#8172b3'
BROWN = '#937860'
PINK = '#da8bc3'
GRAY = '#8c8c8c'
YELLOW = '#ccb974'
CYAN = '#64b5cd'

MIDDLE_BLUE = "#769ad1"
LIGHT_BLUE = '#9bbbea'
LIGHT_GREEN = '#b2df8a'
LIGHT_PURPLE = '#cab2d6'
LIGHT_GRAY = '#b7b7b7'
SUPER_LIGHT_GRAY = '#e3e3e3'
MIDDLE_GRAY = '#a0a0a0'
LIGHT_BROWN = '#e2cbb8ff'

GREATER_EQUAL = u'\u2265'

BCG_PAL = [GREEN, ORANGE, PURPLE, RED, BLUE, PINK, CYAN, BROWN, YELLOW, GRAY,
           LIGHT_PURPLE, LIGHT_GREEN, LIGHT_BLUE, LIGHT_GRAY, MIDDLE_BLUE]

CA = 'Ca'
SA = 'Sa'
LPS = 'LPS'
MT = 'Mt'

NEG = 'Negative'
POS = 'Positive'

VISITS_TO_TIMEPOINTS = {'V1': 'd0', 'V2': 'd14', 'V3': 'd90'}

NAMES = {
    'CYTO': 'cytokine',
    'CM': 'inflaMark',
    'WB_PER_ML': 'cellFreq',
    'PBMC_PERC': 'pbmc%',
    'ATAC': 'ATACseq',
    'SNP': 'SNP',
}

LONG_NAMES = {
    'CYTO': 'Cytokine',
    'CM': 'Inflammatory marker',
    'WB_PER_ML': 'Cell type',
    'PBMC_PERC': 'PBMC cell type',
}

HUMAN_NAMES = {
    'CYTO': 'cytokine production capacity',
    'CM': 'inflammatory marker concentrations',
    'WB_PER_ML': 'immune cell frequencies',
    'PBMC_PERC': 'PBMC cell-type percentages',
}

IMMUNITY_TYPES = {
    'innate_MTB_24h_wo_LAC': 'innate after Mt stimulus',
    'heterologous_nonspecific_7d': 'heterologous immunity',
    'adaptive_MTB_7d': 'adaptive immunity',
    'innate_nonspecific_24h_wo_LAC_IL10_IL1ra': 'trained immunity',
    'innate_nonspecific_24h_wo_LAC': 'trained immunity',
}

SHORT_IMMUNITY_TYPES = {
    'innate_MTB_24h_wo_LAC': 'MTB24H',
    'heterologous_nonspecific_7d': 'HETERO',
    'adaptive_MTB_7d': 'ADAPT',
    'innate_nonspecific_24h_wo_LAC_IL10_IL1ra': 'TRIM',
    'innate_nonspecific_24h_wo_LAC': 'TRIM',
}

EWAS_PHENOTYPES = ['innate_nonspecific_24h_wo_LAC', 'adaptive_MTB_7d',
                   'heterologous_nonspecific_7d', 'innate_MTB_24h_wo_LAC']

RNA_D90 = 'Unstim. day 90 vs. unstim. day 0'
RNA_LPS = 'LPS day 0 vs. unstim. day 0'
RNA_TRIM = 'LPS day 90 vs. LPS day 0'
RNA_TRIM_ABOVE_BCG = 'LPS day 90 vs. unstim. day 90'

PADJ_COL = 'Adjusted P-value'
PVAL_COL = 'P-value'
R_COL = 'Spearman R'
COEF_COL = 'Coefficient'
L2FC_COL = 'Log2 fold change'
T_COL = 'T-value'
F_COL = 'F-value'
ODDS_RATIO_COL = 'Odds ratio'
OVERLAP_COL = 'Overlap'
N_TESTS_COL = 'N tests'
LIBRARY_COL = 'Library'
REGION_SET_COL = 'Region set'
GENE_SET_COL = 'Gene set'
ASSOC_COL = 'Association'
COMPARISON_COL = 'Comparison'
PHENO_COL = 'Phenotype'
REGION_COL = 'Region'
GENE_COL = 'Gene'
NORM_DISP_COL = 'Dispersion (normalized)'
FVE_COL = 'FVE'
CELL_TYPE_COL = 'Cell type'
HOST_FACTOR_COL = 'Host factor'
SCAR_COL = 'Scar size'
GENOMIC_PROPERTY_COL = 'Genomic property'
GENOMIC_LOCATION_COL = 'Genomic location'
FEATURE_TYPE_COL = 'Genomic feature'
REGULATORY_ANNOT_COL = 'Ensembl Regulatory Build'
PROMOTER_MAP_COL = 'Promoter mapping'
DISTAL_MAP_COL = 'Distal mapping'
CLUSTER_COL = 'Cluster'

TSS_PROXIMAL_LABEL = 'TSS proximal'
GENE_BODY_LABEL = 'Gene body'
DISTAL_LABEL = 'Distal (< 10kb)'
INTERGENIC_LABEL = 'Intergenic ({} 10kb)'.format(GREATER_EQUAL)

PROTEIN_CODING = 'Protein-coding'
LNCRNA = 'LncRNA'
IGTR = 'Immunoglobulin (IG) and T-cell receptor (TR)'

CTCF_BINDING = 'CTCF binding site'
TF_BINDING = 'TF binding site'
ENHANCER = 'Distal enhancer'
OPEN_CHROMATIN = 'Open chromatin'
PROMOTER = 'Promoter with TSS'
PROMOTER_FLANK = 'Promoter-flanking region'
UNASSIGNED = 'Unassigned'

RENAME_REG_BUILD = {
    'CTCF_binding_site': CTCF_BINDING,
    'TF_binding_site': TF_BINDING,
    'enhancer': ENHANCER,
    'open_chromatin_region': OPEN_CHROMATIN,
    'promoter': PROMOTER,
    'promoter_flanking_region': PROMOTER_FLANK,
    'reg_NONE': UNASSIGNED
}

KEGG = 'KEGG_2019_Human'
WIKI_PATHWAYS = 'WikiPathways_2019_Human'
BCG_PATHWAYS = 'BCG_pathways'
LM22 = 'LM22'
CHEA = 'ChEA_2022'
GO_BIO_PROCESS = 'GO_Biological_Process_2018'
IMMUNE_SIG_DB = 'ImmuneSigDB_v7.1'
JASPAR = 'JASPAR_2022'
HOCOMOCO = 'HOCOMOCO_v11'
EPI_ROADMAP = 'Roadmap_Epigenomics_r9'

LIB_TO_NAME = {
    KEGG: 'KEGG',
    GO_BIO_PROCESS: 'GO Biological Process',
    JASPAR: 'JASPAR transcription factor binding site motifs',
    HOCOMOCO: 'HOCOMOCO transcription factor binding site motifs',
    EPI_ROADMAP: 'Roadmap Epigenomics'
}

def get_ENR_suppl_cols():
    return ['Term', 'Odds Ratio', 'P-value', 'Adjusted P-value', 'Overlap']


def rename_ENR_suppl_cols(region_sets=False):
    return {'Term': REGION_SET_COL if region_sets else GENE_SET_COL, 'Odds Ratio': ODDS_RATIO_COL,
            'P-value': PVAL_COL, 'Adjusted P-value': PADJ_COL, 'Overlap': OVERLAP_COL}


def get_LR_suppl_cols():
    return ['Coef', 'stat', 'p.value', 'padj']


def rename_LR_suppl_cols(F_test=False, l2fc=False):
    return {'Coef': L2FC_COL if l2fc else COEF_COL, 'stat': F_COL if F_test else T_COL, 'p.value': PVAL_COL, 'padj': PADJ_COL}


def get_limma_suppl_cols(coef, F_test=False, include_coef=True, include_stat=True, include_pval=True, include_padj=True):
    mask = [include_coef, include_stat, include_pval, include_padj]
    return np.asarray(['Coef.{}'.format(coef),
                       '{}.{}'.format('F' if F_test else 't', coef),
                       'p.value.{}'.format(coef),
                       'padj.{}'.format(coef)])[mask].tolist()


def rename_limma_suppl_cols(coef, l2fc=False):
    return {'Coef.{}'.format(coef): L2FC_COL if l2fc else COEF_COL, 'F.{}'.format(coef): F_COL, 't.{}'.format(coef): T_COL,
            'p.value.{}'.format(coef): PVAL_COL, 'padj.{}'.format(coef): PADJ_COL}


#################
# Suppl. Tables #
#################

# S1
COHORT_AND_STATS = 'S01_sample_annotation' # cohort and ATAC-seq statistics
FULL_REGION_ANNOT = 'Ext_S01.1_region_annotation'
# S2
IMMUNE_ANALYSIS = 'S02_immune_functions' # all analysis related to immune function
# S3
ATAC_V1 = 'S03_ATACseq_d0'  # ATAC-seq d0 with variance, host factors, inflam. markers, cytokines, changes in adaptive/trained cytokines
# ExtS3 -- full (all regions) on the suppl. website
ATAC_V1_FULL_HOST_FACTORS = 'Ext_S03.1_ATACseq_d0_hostFactors'  # all regions for host factors
# S4
ENRICH_V1 = 'S04_ATACseq_d0_geneSets'  # ATAC-seq d0 with variance, host factors, inflam. markers, cytokines, changes in adaptive/trained cytokines
# ExtS4 -- promoters and regions sets on the suppl. website
ENRICH_V1_PROMOTERS = 'Ext_S04.1_ATACseq_d0_geneSets_promoters'  # S4 for promoters
ENRICH_V1_LOLA = 'Ext_S04.2_ATACseq_d0_regionSets'
# S8
ATAC_V2_V3 = 'S08_ATACseq_d14FC_d90FC'  # changes in ATAC-seq related to season, trained and adaptive immunity
# ExtS8 -- full (all regions) on the suppl. website
ATAC_V2_V3_FULL = 'Ext_S08.1_ATACseq_d14FC_d90FC'  # all regions for changes related to seasonality
# S9
ENRICH_V2_V3 = 'S09_ATACseq_d14FC_d90FC_geneSets'  # changes in ATAC-seq related to season, trained and adaptive immunity
# ExtS9 -- promoters and regions sets on the suppl. website
ENRICH_V2_V3_PROMOTERS = 'Ext_S09.1_ATACseq_d14FC_d90FC_geneSets_promoters'  # S8 for promoters
# S10
ENRICH_V2_V3_LOLA = 'S10_ATACseq_d14FC_d90FC_regionSets'

DIRECTIONS = {
    -1: NEG,
    1: POS
}
STIMULUS_COLORS = {CA: PURPLE, LPS: PINK, SA: ORANGE, MT: GREEN}
DURATION_COLORS = {'24h': YELLOW, '7d': CYAN}

FONT_SIZE = 6
SMALL_FONT = 4.5
ASTERISK_FONT = 4
DPI = 1000
RASTER = True

RC_PAPER = {
    "font.sans-serif": ["Helvetica"],
    "legend.title_fontsize": SMALL_FONT,
    "figure.titlesize": FONT_SIZE,

    "font.size": FONT_SIZE,
    "axes.labelsize": FONT_SIZE,
    "axes.titlesize": FONT_SIZE,
    "xtick.labelsize": FONT_SIZE,
    "ytick.labelsize": FONT_SIZE,
    "legend.fontsize": FONT_SIZE,

    "axes.linewidth": 0.5,
    "grid.linewidth": 0.5,
    "lines.linewidth": 1,
    "lines.markersize": 3,
    "patch.linewidth": 0.5,

    "xtick.major.width": 0.5,
    "ytick.major.width": 0.5,
    "xtick.minor.width": 0.5,
    "ytick.minor.width": 0.5,

    "xtick.major.size": 2,
    "ytick.major.size": 2,
    "xtick.minor.size": 1,
    "ytick.minor.size": 1,

    "xtick.major.pad": 2,
    "ytick.major.pad": 2,
    "xtick.minor.pad": 2,
    "ytick.minor.pad": 2,
    "axes.labelpad": 2
}
