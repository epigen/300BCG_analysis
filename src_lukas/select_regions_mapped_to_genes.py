import os
from misc import peaks_to_genes, get_peak_annot
from misc import PEAKS_TO_GENES, TSS_PROXIMAL, GENE_AND_DISTAL_10kb, DISTAL_1Mb, METADATA, DATA
from misc import PEAKS_TO_LNC_RNA, LNC_RNA_TSS_PROXIMAL, LNC_RNA_GENE_AND_DISTAL_10kb, LNC_RNA_DISTAL_1Mb


def main():
    df = get_peak_annot(fn=os.path.join(DATA, 'DE', 'peaks_filtered_PBMC.csv.gz'))
    print('all regions', df.shape)
    for region_filter in [TSS_PROXIMAL, GENE_AND_DISTAL_10kb]:
        filtered_df = df.loc[peaks_to_genes(df, **((PEAKS_TO_GENES if region_filter in PEAKS_TO_GENES else PEAKS_TO_LNC_RNA)[region_filter]))].copy()
        print(region_filter, filtered_df.shape)
        filtered_df.to_csv(
            os.path.join(METADATA, 'gene_set_libraries', 'peaks_filtered_PBMC.{}.csv.gz'.format(region_filter)))


if __name__ == '__main__':
    main()
