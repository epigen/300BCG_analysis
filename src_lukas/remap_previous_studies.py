import os
import sys

import pandas as pd
import pybedtools
import misc
import bcg_utils as utils
from misc import DATA
import traceback


def remap_study(study, input_dir, ouput_dir,
                ens_to_gene, ens_to_synonyms, ens_to_gene_mm, ens_to_synonyms_mm,
                ens_mm_to_hs, reference_regions, save_peak_annotation=False, verbose=1, n_threads=1):

    df = pd.read_csv(os.path.join(input_dir, f'{study.name}.csv'))
    print(f"\n-------------------------\nINFO: {study.name} with {study['orig_columns']} and {len(df)} rows\n-------------------------")

    old_index = df.index.copy()
    assert '300BCG_region' not in df.columns or df['300BCG_region'].isnull().all()
    df = df.drop('300BCG_region', axis=1)

    for orig_col in study['orig_columns'].split(','):
        assert not df[orig_col].isnull().any()

    if 'gene' in study['orig_columns'] and study['assay'] not in ['pathway_databases', 'various', 'RNA-seq', 'scRNA-seq']:
        print(f"WARNING: {study.name} assay: {study['assay']} and orig_columns: {study['orig_columns']}")

    # HUMAN
    if study['orig_columns'] in ['hg19_region', 'hg38_region']:
        if study['orig_columns'] == 'hg19_region':
            hg38_regions, unmapped, n_nulls = utils.liftover(df['hg19_region'], remove_tmp=False)
            df.loc[hg38_regions.index, 'hg38_region'] = hg38_regions
            assert unmapped is None or df['hg38_region'].isnull().sum() == len(unmapped) + n_nulls
            if df['hg38_region'].isnull().sum() != 0:
                print(f"WARNING: {study.name} {df['hg38_region'].isnull().sum()}/{len(df)} unmapped hg19 regions")
        df, annot_df = utils.regions_to_genes(
            df=df, region_col='hg38_region', organism=study['organism'], config_fn=misc.UROPA_CONFIG_FN,
            gtf_fn=misc.GENCODE_GTF_FN[study['organism']], prefix=study.name, n_threads=n_threads, return_annot=True
        )
        if save_peak_annotation:
            annot_df.to_csv(os.path.join(ouput_dir, f"{study.name}.peak_annot.csv"), index=False)

        regions_bed = pybedtools.BedTool.from_dataframe(utils.regions_to_bed(df['hg38_region']))
        peaks_300BCG_bed = pybedtools.BedTool.from_dataframe(reference_regions.reset_index()[['chr', 'start', 'end', 'peak_id']])
        for window_size in [0, 500, 1000, 5000, 10000]:
            matched_df = regions_bed.window(peaks_300BCG_bed, w=window_size).to_dataframe()
            matched_df.columns = ['chr', 'start', 'end', 'region_id',
                                  'chr_300BCG', 'start_300BCG', 'end_300BCG', '300BCG_region']
            matched_df['300BCG_region'] = matched_df.groupby('region_id')['300BCG_region'].transform(
                lambda x: ';'.join(sorted(x)))
            matched_df = matched_df[['region_id', '300BCG_region']].drop_duplicates()
            matched_df = matched_df.set_index('region_id')
            assert matched_df.index.is_unique
            col_name = f"300BCG_region{f'_{window_size}bp' if window_size != 0 else ''}"
            df.loc[matched_df.index, col_name] = matched_df['300BCG_region']

    else:
        if study['orig_columns'] == 'human_gene_ensembl':
            assert study['organism'] in ['hs', 'various']
            df = utils.fill_in_gene_names(df, ens_to_gene)
            df = utils.fix_ensembl_ids(df, ens_to_gene, ens_to_synonyms, organism='human', raise_mismatch=True, verbose=verbose, warning_prefix=f'{study.name}: ')

        elif study['orig_columns'] == 'human_gene_name':
            assert study['organism'] in ['hs', 'various']
            assert df['human_gene_ensembl'].isnull().sum() < 0.1 * len(df)
            df = utils.fill_in_ensembl_ids(df, ens_to_gene, organism='human', case_sensitive=False)
            df = utils.fix_ensembl_ids(df, ens_to_gene, ens_to_synonyms, organism='human', verbose=verbose, warning_prefix=f'{study.name}: ')

        elif study['orig_columns'] == 'human_gene_name,human_gene_ensembl':
            assert study['organism'] in ['hs', 'various']
            df = utils.fix_ensembl_ids(df, ens_to_gene, ens_to_synonyms, organism='human', verbose=verbose, warning_prefix=f'{study.name}: ')

        # MOUSE
        elif study['orig_columns'] == 'mm10_region':
            assert study['organism'] == 'mm'
            df, annot_df = utils.regions_to_genes(
                df=df, region_col='mm10_region', organism=study['organism'], config_fn=misc.UROPA_CONFIG_FN,
                gtf_fn=misc.GENCODE_GTF_FN[study['organism']], prefix=study.name, n_threads=n_threads, return_annot=True
            )
            annot_df.to_csv(os.path.join(ouput_dir, f"{study.name}.peak_annot.csv"), index=False)
            df = utils.mouse_to_human(df, ens_mm_to_hs)

        elif study['orig_columns'] == 'mouse_gene_name':
            assert study['organism'] == 'mm'
            assert df['mouse_gene_ensembl'].isnull().sum() < 0.25 * len(df)
            df = utils.fill_in_ensembl_ids(df, ens_to_gene_mm, organism='mouse', case_sensitive=False)
            df = utils.fix_ensembl_ids(df, ens_to_gene_mm, ens_to_synonyms_mm, organism='mouse', verbose=verbose, warning_prefix=f'{study.name}: ')
            df = utils.mouse_to_human(df, ens_mm_to_hs)

        elif study['orig_columns'] in ['mouse_gene_name,mouse_gene_ensembl']:
            assert study['organism'] == 'mm'
            df = utils.fix_ensembl_ids(df, ens_to_gene_mm, ens_to_synonyms_mm, organism='mouse', verbose=verbose, warning_prefix=f'{study.name}: ')
            df = utils.mouse_to_human(df, ens_mm_to_hs)

        else:
            raise ValueError

        df['300BCG_region'] = misc.genes_to_regions(
            gene_ids=df['human_gene_ensembl'],
            peaks_and_genes=reference_regions['gene_id'].str.rsplit('.', n=1, expand=True)[0]
        )

    assert df.index.equals(old_index)

    return df


def main():
    study_names = sys.argv[1:] if len(sys.argv) > 1 else None

    with_synonyms = 'r97' not in misc.HUMAN_ENSEMBL_IDS_FN

    ens_to_gene, ens_to_synonyms = utils.get_ens_to_gene_map(misc.HUMAN_ENSEMBL_IDS_FN, with_synonyms=with_synonyms)
    ens_to_gene_mm, ens_to_synonyms_mm = utils.get_ens_to_gene_map(misc.MOUSE_ENSEMBL_IDS_FN, with_synonyms=with_synonyms)
    ens_mm_to_hs = utils.get_ens_mm_to_hs_map(misc.ENSEMBL_MOUSE_ORTHOLOGS_FN, check=ens_to_gene_mm)
    peaks_300BCG_df = misc.get_peak_annot()

    meta_fn = os.path.join(DATA, 'trained_immunity_gene_sets_20221007', 'metadata_trainedImmunityGenesRegions.csv')
    input_dir = os.path.join(DATA, 'trained_immunity_gene_sets_20221007', 'sets')
    ouput_dir = utils.make_dir(DATA, 'trained_immunity_gene_sets_20221007', 'converted_r97')

    meta_df = pd.read_csv(meta_fn, index_col=0)
    assert not meta_df['PMID'].isnull().any()

    SAVE = True
    SAVE_PEAK_ANNOT = True
    N_THREADS = 16
    EACH_ONCE = False
    tested = []

    for idx, study in meta_df.iterrows():

        if study_names is not None and study.name not in study_names:
            continue

        if EACH_ONCE and study['orig_columns'] in tested:
            continue
        else:
            tested.append(study['orig_columns'])

        try:
            df = remap_study(study, input_dir=input_dir, ouput_dir=ouput_dir,
                             ens_to_gene=ens_to_gene, ens_to_synonyms=ens_to_synonyms,
                             ens_to_gene_mm=ens_to_gene_mm, ens_to_synonyms_mm=ens_to_synonyms_mm,
                             ens_mm_to_hs=ens_mm_to_hs, reference_regions=peaks_300BCG_df,
                             save_peak_annotation=SAVE_PEAK_ANNOT, n_threads=N_THREADS)
            if df is None:
                raise ValueError

            if SAVE:
                df.to_csv(os.path.join(ouput_dir, f"{study.name}.csv"), index=False)
        except Exception as e:
            print(f"ERROR {study.name}: {e}")
            traceback.print_exception(*sys.exc_info())


if __name__ == '__main__':
    main()
