import matplotlib as mpl
mpl.use('Agg')
from misc import *


USE_BINARY = False

suffix = ''
print(suffix)
data_folder = '{}'.format(DATA)
print(data_folder)
analysis_folder = os.path.join(data_folder, "analysis")
sample_annotations_file = os.path.join(analysis_folder, "complete_metadata.csv")
print(sample_annotations_file)
peaks_annotations_file = os.path.join(analysis_folder, "peaks_characterization{suffix}.csv.gz".format(suffix=suffix))
binary_file = os.path.join(analysis_folder, "quantification_binary{suffix}.csv.gz".format(suffix=suffix))
count_file = os.path.join(analysis_folder, "quantification{suffix}.csv.gz".format(suffix=suffix))
print(count_file)

sample_annotations = pd.read_csv(sample_annotations_file, index_col='DEMUX:BIOSAMPLE')
print('sample_annotations', sample_annotations.shape)
sample_annotations = sample_annotations[sample_annotations['QC:PASS']]
print('sample_annotations (QC:PASS)', sample_annotations.shape)

peaks_annotations = pd.read_csv(peaks_annotations_file, index_col='peak_id')
print('\nChromosomes:', sorted(set(peaks_annotations['chr'])))
print('peaks_annotations', peaks_annotations.shape)

chromosomes = peaks_annotations['chr'].isin(['chr{}'.format(i) for i in range(1, 23)])
selected_peaks = peaks_annotations[chromosomes]
print('selected chromosomes:', sorted(set(selected_peaks['chr'])))
print('selected_peaks', selected_peaks.shape)

not_selected = peaks_annotations.loc[~peaks_annotations.index.isin(selected_peaks.index)]
not_selected = not_selected.loc[chromosomes]
print('\nNot selected:')
print(not_selected.groupby(['characterization', 'feat_type'])['gene_name'].count())

counts = pd.read_csv(count_file, index_col='ID').loc[selected_peaks.index, sample_annotations.index]
if USE_BINARY:
    binary_peaks_df = pd.read_csv(binary_file, index_col='ID').loc[selected_peaks.index, sample_annotations.index]
cpm = tpm(counts, transcript_len=1, norm_factors=1, log=False, pseudocount=0.5, libsize_psuedocount=0, samples_as_rows=False)

assert counts.columns.equals(sample_annotations.index)
assert cpm.columns.equals(sample_annotations.index)
assert cpm.columns.equals(counts.columns)
assert cpm.index.equals(counts.index)

min_count_mask = {}
for celltype, visits in [
    ('PBMC', ['V1', 'V2', 'V3']),
    # ('nkcell', ['V1', 'V2', 'V3']),
    # ('monocyte', ['V1', 'V2', 'V3']),
    # ('cd8t', ['V1', 'V2', 'V3'])
]:
    min_count_mask[celltype] = {}
    for large_min_n_proportion in [0.1 if USE_BINARY else 0.7]:
        min_group_ns = pd.get_dummies(sample_annotations.loc[(sample_annotations['SAMPLE:TISSUE'] == celltype) & (sample_annotations['SAMPLE:VISIT'].isin(visits)), 'SAMPLE:VISIT']).sum()
        print('\n', celltype, 'min_group_ns\n', min_group_ns)

        if USE_BINARY:
            min_count_mask[celltype][large_min_n_proportion] = \
                (binary_peaks_df.loc[:, sample_annotations['SAMPLE:TISSUE'] == celltype] != 0).sum(axis=1) >= (min_group_ns.min() * large_min_n_proportion)
        else:
            min_count_mask[celltype][large_min_n_proportion] = filter_by_reads(
                counts.loc[:, (sample_annotations['SAMPLE:TISSUE'] == celltype) & (sample_annotations['SAMPLE:VISIT'].isin(visits))],
                min_group_n=min_group_ns.min(),
                cpm_df=cpm.loc[:, (sample_annotations['SAMPLE:TISSUE'] == celltype) & (sample_annotations['SAMPLE:VISIT'].isin(visits))],
                large_min_n_proportion=large_min_n_proportion,
                for_large_min_n_use_only_proportion=True,
                samples_as_rows=False)

        print('before filtering',
              counts.loc[:, (sample_annotations['SAMPLE:TISSUE'] == celltype) & (sample_annotations['SAMPLE:VISIT'].isin(visits))].shape)
        print('after filtering',
              counts.loc[min_count_mask[celltype][large_min_n_proportion], (sample_annotations['SAMPLE:TISSUE'] == celltype) & (sample_annotations['SAMPLE:VISIT'].isin(visits))].shape)

if USE_BINARY:
    del binary_peaks_df

cpm = np.log2(cpm)

debug = False
for cell_types, visits in [
    (['PBMC'], ['V1', 'V2', 'V3']),
    (['cd8t', 'monocyte', 'nkcell'], ['V1']),
    # (['cd8t'], ['V1', 'V2', 'V3']),
    # (['monocyte'], ['V1', 'V2', 'V3']),
    # (['nkcell'], ['V1', 'V2', 'V3'])
]:
    cell_annotations = sample_annotations.loc[sample_annotations['SAMPLE:TISSUE'].isin(cell_types) & sample_annotations['SAMPLE:VISIT'].isin(visits)]
    cell_types = '_'.join(sorted(cell_types))
    if '_' in cell_types:
        _min_count_mask = min_count_mask['PBMC']
    else:
        _min_count_mask = min_count_mask[cell_types]
    print(cell_types, cell_annotations.shape[0])
    for large_min_n_proportion in _min_count_mask:
        print('Saving the data for', cell_types, 'with large_min_n_proportion of', large_min_n_proportion)
        de_folder = make_dir(
            data_folder, ('DE{suffix}_new' if len(_min_count_mask) == 1 else 'DE{}_{}'.format('{suffix}', large_min_n_proportion)).format(suffix='{}{}'.format(suffix, '_binary' if USE_BINARY else '')))

        pd.Series(cpm.loc[_min_count_mask[large_min_n_proportion]].index).to_csv(os.path.join(de_folder, 'peaks_{}.csv'.format(cell_types)), index=False, header=False)

        fig, axs = plt.subplots(2, 2, figsize=(10, 10))
        plt.subplots_adjust(wspace=0.5, hspace=0.5)

        ax = sns.scatterplot(cpm.loc[:, cell_annotations.index].mean(axis=1),
                             cpm.loc[:, cell_annotations.index].std(axis=1),
                             ax=axs[0, 0], rasterized=True)
        ax.set_xlabel('Mean')
        ax.set_ylabel('Std')
        sns.despine()
        ax.set_title('{} (before filtering)\n{} peaks'.format(cell_types, cpm.shape[0]))

        ax = sns.scatterplot(cpm.loc[_min_count_mask[large_min_n_proportion], cell_annotations.index].mean(axis=1),
                             cpm.loc[_min_count_mask[large_min_n_proportion], cell_annotations.index].std(axis=1),
                             ax=axs[0, 1], rasterized=True)
        ax.set_xlabel('Mean')
        ax.set_ylabel('Std')
        sns.despine()
        ax.set_title('{} (after filtering)\n{} peaks'.format(cell_types, _min_count_mask[large_min_n_proportion].sum()))

        plot_sequenced_samples(cpm.loc[:, cell_annotations.index], n_samples=None, ax=axs[1, 0], xlabel='Log$2$ CPM', ylabel='Density', title='{} (before filtering)'.format(cell_types), samples_as_rows=False)
        plot_sequenced_samples(cpm.loc[_min_count_mask[large_min_n_proportion], cell_annotations.index], n_samples=None, ax=axs[1, 1], xlabel='Log$2$ CPM', ylabel='Density', title='{} (after filtering)'.format(cell_types), samples_as_rows=False)
        
        savefig(os.path.join(de_folder, 'filtering_{}.pdf'.format(cell_types)), dpi=300)
        plt.close()

        if not debug:
            _selected_counts = counts.loc[_min_count_mask[large_min_n_proportion], cell_annotations.index]
            _norm_factors = calc_norm_factors(_selected_counts, method='TMM')
            _norm_factors = pd.Series(_norm_factors, index=_selected_counts.columns, name='norm_factor')
            _selected_cpm = tpm(_selected_counts, norm_factors=_norm_factors, log=True, pseudocount=0.5,
                                 libsize_psuedocount=1)

            peaks_annotations.loc[_selected_counts.index].to_csv(
                os.path.join(de_folder, "peaks_filtered_{}.csv.gz".format(cell_types)), index_label='peak_id')

            assert _selected_counts.columns.str.endswith('_ATAC_R1').all()
            _selected_counts.columns = _selected_counts.columns.str.replace('_ATAC_R1$', '')
            _selected_counts.to_csv(os.path.join(de_folder, "quantification_filtered_{}.csv.gz".format(cell_types)),
                                    index_label='ID')

            _cell_annotations = cell_annotations.copy()
            assert _cell_annotations.index.str.endswith('_ATAC_R1').all()
            _cell_annotations.index = _cell_annotations.index.str.replace('_ATAC_R1$', '')
            _cell_annotations.to_csv(os.path.join(de_folder, "samples_filtered_{}.csv".format(cell_types)),
                                     index_label='DEMUX:BIOSAMPLE')

            assert _selected_cpm.columns.str.endswith('_ATAC_R1').all()
            _selected_cpm.columns = _selected_cpm.columns.str.replace('_ATAC_R1$', '')
            _selected_cpm.to_csv(os.path.join(de_folder, "normalized_log2_CPM_{}.csv.gz".format(cell_types)),
                                  index_label='ID')

            assert _norm_factors.index.str.endswith('_ATAC_R1').all()
            _norm_factors.index = _norm_factors.index.str.replace('_ATAC_R1$', '')
            _norm_factors.to_csv(os.path.join(de_folder, "norm_factors_{}.csv.gz".format(cell_types)), index_label='ID')

    if not debug:
        pbmc = pd.read_csv(os.path.join(de_folder, "normalized_log2_CPM_PBMC.csv.gz"), index_col=0)
        others = pd.read_csv(os.path.join(de_folder, "normalized_log2_CPM_cd8t_monocyte_nkcell.csv.gz"), index_col=0)
        data = pd.concat([pbmc, others], axis=1).T
        pca_data = PCA(n_components=2).fit_transform(StandardScaler().fit_transform(data))
        sns.scatterplot(pca_data[:, 0], pca_data[:, 1], hue=data.index.str.split('_', expand=True).get_level_values(2))
        savefig(os.path.join(de_folder, 'PCA_scatterplot.pdf'), dpi=300)
        plt.close()
