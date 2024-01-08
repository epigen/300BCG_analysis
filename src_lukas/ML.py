#!/usr/bin/env python
# coding: utf-8

import matplotlib as mpl
import matplotlib.pyplot as plt
import joblib
import gzip
mpl.use('Agg')
from sklearn.cluster import KMeans
from misc import *
import bcg_utils

from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
from umap import UMAP

sns.set_color_codes()
sns.set_context('talk')
sns.set_style('ticks')
mpl.rcParams['legend.title_fontsize'] = 16
mpl.rcParams['figure.titlesize'] = 18

# assert os.environ['PYTHONHASHSEED'] == '0'
os.environ["PYTHONWARNINGS"] = "ignore:The max_iter was reached:UserWarning"


def transform_features(X_df, config, save_data=False):

    X_df = X_df.loc[~X_df.isnull().any(axis=1)]

    if 'min_mean' in config and config['min_mean'] is not None:
        X_df = X_df.loc[:, X_df.mean(axis=0) >= config['min_mean']]
        print('min_mean X', X_df.shape)

    if 'max_mean' in config and config['max_mean'] is not None:
        X_df = X_df.loc[:, X_df.mean(axis=0) <= config['max_mean']]
        print('max_mean X', X_df.shape)

    if 'highly_variable' in config and config['highly_variable'] is not None:
        import scanpy as sc
        from anndata import AnnData
        adata = AnnData(X=np.log(np.power(2, X_df.astype(float)) + 1))
        hvr_df = sc.pp.highly_variable_genes(adata, n_top_genes=config['highly_variable'], flavor='seurat',
                                             inplace=False)
        hvr_df.index = adata.var_names
        print(fn)
        sc.pl.highly_variable_genes(hvr_df, show=False, save='high_var_regions.pdf')
        X_df = X_df.loc[:, hvr_df['highly_variable']]
        del adata
        del hvr_df
        print('highly_variable X', X_df.shape)

    if config.get('global_PCA') is not None or config.get('global_KPCA') is not None or config.get(
            'global_UMAP') is not None:
        if (config['X_name'] in ['CYTO', 'CM']):
            X_df = pd.DataFrame(data=IterativeImputer(max_iter=10, random_state=config['seed']).fit_transform(X_df),
                                index=X_df.index, columns=X_df.columns)

        _PCA_COLS = r'^CONS0|^CYTO:|^CM:|^PBMC_PER_ML:'
        if config.get('global_PCA') is not None:
            _feature_prefix = 'PC'
            print('PCA')
            _pca = PCA(n_components=config['global_PCA'] if config['global_PCA'] > 0 else None,
                            random_state=config['seed'], svd_solver='full')
            _columns = X_df.columns[X_df.columns.str.contains(_PCA_COLS)]
            embedding = _pca.fit_transform(
                StandardScaler().fit_transform(X_df.loc[:, _columns]))
            if save_data:
                loadings = pd.DataFrame(data=_pca.components_,
                                        index=['{}{}'.format(_feature_prefix, i + 1) for i in range(_pca.components_.shape[0])],
                                        columns=_columns)
                loadings.to_csv(os.path.join(ML_RESULTS_ROOT, 'PCA_loadings.csv.gz'))
        elif config.get('global_KPCA') is not None:
            embedding = KernelPCA(kernel='poly', degree=config['global_KPCA_degree'], n_components=config['global_KPCA'],
                                  random_state=config['seed'], eigen_solver='dense').fit_transform(
                StandardScaler().fit_transform(X_df.loc[:, X_df.columns.str.contains(_PCA_COLS)]))
            _feature_prefix = 'KPC'
        elif config.get('global_UMAP') is not None:
            embedding = UMAP(n_neighbors=config['n_neighbors'], min_dist=config['min_dist'],
                             n_components=config['global_UMAP'], random_state=config['seed']).fit_transform(
                StandardScaler().fit_transform(X_df.loc[:, X_df.columns.str.contains(_PCA_COLS)]))
            _feature_prefix = 'UMAP'

        embedding = pd.DataFrame(data=embedding, index=X_df.index,
                                 columns=['{}{}'.format(_feature_prefix, i + 1) for i in range(embedding.shape[1])])
        X_df = pd.concat([X_df.loc[:, ~X_df.columns.str.contains(_PCA_COLS)], embedding], axis=1)
        print('X', X_df.shape, X_df.columns)

    if config.get('poly') is not None:
        _POLY_COLS = r'^PC|^KPC|^UMAP|^CONS0|^CYTO:|^CM:|^PBMC_PER_ML:'
        from sklearn.preprocessing import PolynomialFeatures
        poly = PolynomialFeatures(degree=config['poly'], interaction_only=False, include_bias=False)
        poly_features = poly.fit_transform(X_df.loc[:, X_df.columns.str.contains(_POLY_COLS)])
        poly_features = pd.DataFrame(data=poly_features, index=X_df.index, columns=poly.get_feature_names(X_df.columns[X_df.columns.str.contains(_POLY_COLS)]))
        X_df = pd.concat([X_df.loc[:, ~X_df.columns.str.contains(_POLY_COLS)], poly_features], axis=1)
        print('Polynomial features')
        print(X_df.columns.tolist())
        print('X', X_df.shape, X_df.columns)

    if config['scale']:
        print('Global scaling')
        X_df = pd.DataFrame(data=StandardScaler().fit_transform(X_df), index=X_df.index, columns=X_df.columns)
        # _season_cols = X_df.columns.str.contains(r'^DONOR:IC_DATE_2PI_SIN$|^DONOR:IC_DATE_2PI_COS$')
        # if _season_cols.all():
        #     X_df = pd.DataFrame(data=X_df / np.sqrt(0.5), index=X_df.index, columns=X_df.columns)
        # else:
        #     print('Season and other features')
        #     X_df = pd.concat([pd.DataFrame(data=X_df.loc[:, _season_cols] / np.sqrt(0.5), index=X_df.index, columns=X_df.columns[_season_cols]),
        #                       pd.DataFrame(data=StandardScaler().fit_transform(X_df.loc[:, ~_season_cols]), index=X_df.index, columns=X_df.columns[~_season_cols])],
        #                      axis=1)

    if save_data:
        X_df.to_csv(os.path.join(ML_RESULTS_ROOT, 'X_df.csv.gz'))
    return X_df


def _nested_cv(nested_FS, estimator, X, y, y_real, X_extra, y_extra, groups, samples, peaks, target, test_kfold, cv_kfold,
               param_grid, scoring, alt_scoring, permute_labels, subsample, scale, binarized, merged_folds_scoring, feature_selection_name, feature_selection_mask,
               feature_selection_score, n_features, random_features, n_pcs, pca_first, n_jobs, random_state, verbose, out):

    nested_cv_params = dict(estimator=estimator, X=X, y=y, groups=groups,
                            samples=samples, target=target, test_kfold=test_kfold, cv_kfold=cv_kfold,
                            param_grid=param_grid, scoring=scoring, alt_scoring=alt_scoring,
                            permute_labels=permute_labels, subsample=subsample, scale=scale,
                            feature_selection_score=feature_selection_score, n_features=n_features,
                            n_pcs=n_pcs, n_jobs=n_jobs, random_state=random_state, verbose=verbose, out=out)

    use_regression_and_extras = feature_selection_name is not None and 'regression' in feature_selection_name
    print('use_regression_and_extras:', use_regression_and_extras)

    if nested_FS:
        #assert feature_selection_name is not None, 'Maybe you should switch OFF nested_FS?'
        if feature_selection_name and ('blind_DE' in feature_selection_name or 'rand_DE' in feature_selection_name):
            raise ValueError
            print('NESTED FEATURE SELECTION WITH BLIND DE')
            print('random_features:', random_features)
            return bcg_utils.nested_cv_with_blind_DE(**nested_cv_params, random_features=random_features, peaks=peaks)
        else:
            print('CLASSIC NESTED FEATURE SELECTION')
            print('merged_folds_scoring:', merged_folds_scoring)
            return bcg_utils.nested_cv(
                binarized=binarized.values if binarized is not None else None, strictly_use_proba=True, merged_folds_scoring=merged_folds_scoring, pca_first=pca_first,
                **nested_cv_params
            )
    else:
        raise ValueError
        print('FROZEN FEATURE SELECTION')
        return bcg_utils.frozen_nested_cv(**nested_cv_params,
                                y_real=y_real, X_extra=X_extra, y_extra=y_extra,
                                feature_selection_name=feature_selection_name,
                                feature_selection_mask=feature_selection_mask,
                                add_extra_samples_for_fs=use_regression_and_extras,
                                use_y_real_for_fs=use_regression_and_extras)


show = False
config_list_fn = sys.argv[1]
n_jobs = int(sys.argv[2]) if len(sys.argv) > 2 else None
save_data = (sys.argv[3] == '--save_data') if len(sys.argv) > 3 else False
task_id = int(os.environ['SLURM_ARRAY_TASK_ID'])
print(config_list_fn, n_jobs, task_id)

with open(os.path.join(ML_RESULTS_ROOT, config_list_fn) if not config_list_fn.startswith(ML_RESULTS_ROOT + os.path.sep) else config_list_fn) as f:
    config_fn = f.readlines()[task_id - 1].strip()
print(config_fn)

with open(config_fn, 'r') as f:
    config = yaml.load(f, Loader=yaml.loader.FullLoader)

results_dir = os.path.join(ML_RESULTS_ROOT, get_data_name_from_config(config))
fn = get_filename_from_config(config, results_dir)
log = open(fn.format(data='results', ext='log'), 'w')
print(fn, file=log)


# X
model = get_model_for_visits(config['model'], config['X_visits'])
print('Model for X:', config['X_visits'], model)
data_X_fn = os.path.join(ML_RESULTS_ROOT, model, '{}.{}.hdf'.format(PROJECT, model))
sample_annot_df = pd.read_hdf(data_X_fn, key='{}/sample_annot'.format(config['celltype']))
sample_annot_df.loc[sample_annot_df['SAMPLE:VISIT'] == 'V1', 'DONOR:SCALED_AGE'] = StandardScaler().fit_transform(
    sample_annot_df.loc[sample_annot_df['SAMPLE:VISIT'] == 'V1', 'DONOR:AGE'].values.reshape(-1, 1))
sample_annot_df.loc[sample_annot_df['SAMPLE:VISIT'] == 'V1', 'DONOR:SCALED_BMI'] = StandardScaler().fit_transform(
    sample_annot_df.loc[sample_annot_df['SAMPLE:VISIT'] == 'V1', 'DONOR:BMI'].values.reshape(-1, 1))
sample_annot_df.loc[sample_annot_df['SAMPLE:VISIT'] == 'V1', 'DONOR:SCALED_IC_TIME_REAL'] = StandardScaler().fit_transform(
    sample_annot_df.loc[sample_annot_df['SAMPLE:VISIT'] == 'V1', 'DONOR:IC_TIME_REAL'].values.reshape(-1, 1))
X_df = pd.DataFrame()
for X in config['X']:
    if X == 'ATAC':
        X_df = pd.concat([X_df,
                          pd.read_hdf(data_X_fn, key='{}/atac'.format(config['celltype']))],
                         axis=1)
    else:
        X_df = pd.concat([X_df,
                          bcg_utils.force_df_to_numeric(
                              sample_annot_df.loc[:, sample_annot_df.columns.str.contains(X)], include_binary=True)],
                         axis=1)

        for to_impute_cols in [r'^CYTO:', r'^CM:']:
            _cyto_cols = X_df.columns.str.contains(to_impute_cols)
            if _cyto_cols.sum() != 0:
                print('Imputing.')
                _all_cyto_null = X_df.loc[:, _cyto_cols].isnull().all(axis=1)
                X_df.loc[~_all_cyto_null, _cyto_cols] = impute_cytokines(
                    X_df.loc[~_all_cyto_null, _cyto_cols], max_iter=10, random_state=config['seed'])

X_df = X_df.astype(float)

if config['X_name'] in ['LOG2_PBMC_PERC', 'LOG2_PBMC_PER_ML', 'LOG2_WB_PER_ML']:
    X_df = X_df.loc[:, X_df.columns.str.contains(':MONO$|:MONO/INTERMEDIATE$|:MONO/NON_CLASSICAL$|:T/CD8$|:T/CD4$|:T/CD4/TREG$|:B$|:NK$|:NKT$|:BASO$|:NEUTRO$')]
    print('Transforming with log2:', X_df.columns.values)
    X_df = X_df - X_df.min().min() + 0.1
    assert not X_df.isnull().any().any()
    assert (X_df > 0).all().all()
    X_df = np.log2(X_df)

print('These are the features:\n{}\n'.format(
    ', '.join((['ATAC'] if X_df.columns.str.startswith('CONS0').any() else [])
              + (['PCA_ATAC'] if X_df.columns.str.startswith('PC').any() else [])
              + X_df.columns[~X_df.columns.str.contains('^CONS0|^PC')].tolist())
))

# Y
model = get_model_for_visits(config['model'], config['Y_visits'])
print('Model for Y:', config['Y_visits'], model)
data_Y_fn = os.path.join(ML_RESULTS_ROOT, model, '{}.{}.hdf'.format(PROJECT, model))
_sample_annot_df_index = sample_annot_df.index.copy()
sample_annot_df = pd.read_hdf(data_Y_fn, key='{}/sample_annot'.format(config['celltype']))
assert _sample_annot_df_index.equals(sample_annot_df.index)

if config['Y'].startswith('ATAC'):
    _, _celltype, _top_n = config['Y'].split('_')
    _top_n = int(_top_n)
    _atac_df = pd.read_hdf(data_Y_fn, key='{}/atac'.format(_celltype))
    de_fn = os.path.join(ML_RESULTS_ROOT, 'DE', 'de_{data}_{celltype}.{model}.csv.gz')
    regions = set([])
    for coef in ['V2vsV1', 'V3vsV1']:
        pval_col, coef_col = 'p.value.{coef}'.format(coef=coef), 'Coef.{coef}'.format(coef=coef)
        df = pd.read_csv(de_fn.format(data='results_p5', celltype=_celltype, model='donor'), usecols=[pval_col, coef_col])
        df[pval_col] = -np.log10(df[pval_col])
        df.loc[df[coef_col] < 0, pval_col] *= -1
        regions = regions.union(df[pval_col].sort_values().iloc[:_top_n].index)
        regions = regions.union(df[pval_col].sort_values().iloc[-_top_n:].index)
    Y_df = _atac_df[sorted(regions)]

    sample_annot_df.index = sample_annot_df.index.str.replace('_{}_'.format(config['celltype']), '_')
    X_df.index = X_df.index.str.replace('_{}_'.format(config['celltype']), '_')
    Y_df.index = Y_df.index.str.replace('_{}_'.format(config['celltype']), '_')

    _idx = sorted(set(X_df.index).intersection(Y_df.index))
    sample_annot_df = sample_annot_df.loc[_idx]
    X_df = X_df.loc[_idx]
    Y_df = Y_df.loc[_idx]
else:
    _specials = get_basic_binary_specials()
    if config['binarize'] and not config['scoring'].startswith('quantile_'):
        _flu_season = lambda x: 0 if 6 <= x < 10 else (1 if x < 4 or x >= 11 else np.nan)
        _specials['SAMPLE:VISIT_MONTH_REAL'] = lambda s: s.apply(_flu_season)
        _specials['DONOR:IC_MONTH_REAL'] = lambda s: s.apply(_flu_season)

    if 'IMPUTED_CYTO' in config['Y']:
        assert '|' not in config['Y']
        _regex = '{}|{}'.format(config['Y'], config['Y'].replace('IMPUTED_CYTO', 'CYTO'))
        _real_targets = get_numeric_targets(sample_annot_df, regex=config['Y'], specials=_specials, include_binary=True)
    else:
        _regex = config['Y']
        _real_targets = None
    Y_df = get_numeric_targets(sample_annot_df, regex=_regex, specials=_specials, include_binary=True)
    evening_donors = set(sample_annot_df.loc[sample_annot_df['DONOR:IC_TIME_REAL'] > 13, 'SAMPLE:DONOR'])

print('X', X_df.shape, 'Y', Y_df.shape, file=log)
assert X_df.index.equals(Y_df.index) and X_df.index.equals(sample_annot_df.index)
sample_annot_df = sample_annot_df.set_index(['SAMPLE:DONOR', 'SAMPLE:VISIT'], drop=False, append=False)
X_df.index = sample_annot_df.index
Y_df.index = sample_annot_df.index

if config.get('bootstrap'):
    print('Bootstrapping')
    X_df = X_df.sample(frac=1, replace=True, random_state=config['seed'])
    Y_df = Y_df.loc[X_df.index]
    sample_annot_df = sample_annot_df.loc[X_df.index]
    assert X_df.index.equals(Y_df.index) and X_df.index.equals(sample_annot_df.index)
    new_idx = pd.MultiIndex.from_arrays([['{}:{}'.format(x, i + 1) for i, x in enumerate(X_df.index.get_level_values('SAMPLE:DONOR'))],
                                         X_df.index.get_level_values('SAMPLE:VISIT')],
                                        names=['SAMPLE:DONOR', 'SAMPLE:VISIT'])
    X_df.index = Y_df.index = sample_annot_df.index = new_idx

if '-' not in config['X_visits'][0]:
    print('Transforming BEFORE selecting visits:', config['X_visits'])
    X_df = transform_features(X_df, config, save_data=save_data)
    if save_data:
        print('Data saved, exiting...')
        sys.exit()

assert len(config['X_visits']) == len(config['Y_visits'])
assert set(config['X_visits']) == set(config['Y_visits']) or len(config['X_visits']) == 1
assert all(['-' not in v for v in config['X_visits']]) or len(config['X_visits']) == 1
assert all(['-' not in v for v in config['Y_visits']]) or len(config['Y_visits']) == 1

print('Selecting visits')
X_df = select_visits(df=X_df, visits=config['X_visits'], no_diff_cols=NO_DIFF_COLS)
Y_df = select_visits(df=Y_df, visits=config['Y_visits'], no_diff_cols=NO_DIFF_COLS)
if len(config['X_visits']) == 1:
    sample_annot_df = sample_annot_df.loc[
        sample_annot_df['SAMPLE:VISIT'] == ('V1' if '-' in config['X_visits'][0] else config['X_visits'][0])]
    sample_annot_df.index = sample_annot_df.index.get_level_values('SAMPLE:DONOR')
print('X', X_df.shape, 'Y', Y_df.shape)

print('Removing nulls')
X_df = X_df.loc[~X_df.isnull().any(axis=1)]
Y_df = Y_df.loc[~Y_df.isnull().any(axis=1)]
print('X', X_df.shape, 'Y', Y_df.shape)

# fix indexes
if '-' in config['X_visits'][0] and '-' in config['Y_visits'][0]:
    pass
elif '-' in config['X_visits'][0]:
    assert X_df.index.name == 'SAMPLE:DONOR'
    assert Y_df.index.names == ['SAMPLE:DONOR', 'SAMPLE:VISIT']
    assert set(Y_df.index.get_level_values('SAMPLE:VISIT')) == set(config['Y_visits'])
    Y_df.index = Y_df.index.get_level_values('SAMPLE:DONOR')
elif '-' in config['Y_visits'][0]:
    assert Y_df.index.name == 'SAMPLE:DONOR'
    assert X_df.index.names == ['SAMPLE:DONOR', 'SAMPLE:VISIT']
    assert set(X_df.index.get_level_values('SAMPLE:VISIT')) == set(config['X_visits'])
    X_df.index = X_df.index.get_level_values('SAMPLE:DONOR')
else:
    assert X_df.index.names == ['SAMPLE:DONOR', 'SAMPLE:VISIT']
    assert Y_df.index.names == ['SAMPLE:DONOR', 'SAMPLE:VISIT']
    if len(config['X_visits']) == 1:
        assert set(X_df.index.get_level_values('SAMPLE:VISIT')) == set(config['X_visits'])
        assert set(Y_df.index.get_level_values('SAMPLE:VISIT')) == set(config['Y_visits'])
        X_df.index = X_df.index.get_level_values('SAMPLE:DONOR')
        Y_df.index = Y_df.index.get_level_values('SAMPLE:DONOR')
assert X_df.index.is_unique and Y_df.index.is_unique and sample_annot_df.index.is_unique

X_df, Y_df = X_df.align(Y_df, join='inner', axis=0)
sample_annot_df = sample_annot_df.loc[X_df.index]
print('X', X_df.shape, 'Y', Y_df.shape)
print('X', X_df.shape, 'Y', Y_df.shape, file=log)
assert X_df.index.equals(Y_df.index) and X_df.index.equals(sample_annot_df.index)

if len(config['X_visits']) > 1:
    donors = select_visits(df=sample_annot_df['SAMPLE:DONOR'], visits=config['X_visits'])
elif 'group_by_date' in config and config['group_by_date'] is not None:
    print('group_by_date:', config['group_by_date'])
    if config['group_by_date'] == -1:
        donors = sample_annot_df['DONOR:IC_DATE'].copy()
        print('Grouping based on the actual IC_DATE:', len(set(donors)))
    else:
        clusters = pd.Series(
            ['C{}'.format(c) for c in KMeans(n_clusters=config['group_by_date'], random_state=RANDOM_STATE).fit_predict(sample_annot_df[['DONOR:IC_DATE_REAL']])],
            index=sample_annot_df.index, name='vaccination_date_clusters')
        # sns.scatterplot(x=sample_annot_df['DONOR:IC_DATE_REAL'],
        #                 y=sample_annot_df['thm.innate_nonspecific_24h_wo_LAC_IL10_IL1ra_V3'],
        #                 hue=clusters).legend(bbox_to_anchor=(1, 1))
        # sns.despine()
        # bcg_utils.savefig(fn.format(data='group_by_date', ext='pdf'))
        # plt.show() if show else plt.close()
        donors = clusters
else:
    donors = None
assert donors is None or X_df.index.equals(donors.index)

if '-' in config['X_visits'][0]:
    print('Transforming AFTER selecting visits:', config['X_visits'])
    X_df = transform_features(X_df, config, save_data=save_data)
    if save_data:
        print('Data saved, exiting...')
        sys.exit()

if config['features'] is not None or config['pathways'] is not None or config['DE_filter'] is not None:
    assert config['X'] == ['ATAC']
    peak_annot_df = pd.read_hdf(data_X_fn, key='{}/peak_annot'.format(config['celltype']))
    assert peak_annot_df.index.equals(X_df.columns)

if config['features'] is not None:
    features_mask = peak_annot_df[config['features']['type']] == config['features']['value']
    X_df = X_df.loc[:, features_mask]
    peak_annot_df = peak_annot_df.loc[features_mask]
    print(config['features']['type'], config['features']['value'], 'X', X_df.shape, 'Y', Y_df.shape, file=log)

if config['pathways'] is not None:
    library = bcg_utils.gene_set_library(os.path.join(METADATA, 'gene_set_libraries', config['pathways']['type'] + '.gmt'))
    if config['pathways']['value'] in [config['pathways']['type'], 'all_pathways']:
        genes = [g for p in library for g in library[p]]
        print('Using all pathways in', config['pathways']['type'], ':', len(genes))
    else:
        genes = library[config['pathways']['value']]
        print('Using only', config['pathways']['value'], 'in', config['pathways']['type'], ':', len(genes))
    if config['pathways']['type'].endswith('_regions'):
        assert all([g.startswith('CONS') for g in genes])
        pathways_mask = peak_annot_df.index.isin(genes)
    else:
        pathways_mask = peak_annot_df['gene_name'].isin(genes).values
    X_df = X_df.loc[:, pathways_mask]
    peak_annot_df = peak_annot_df.loc[pathways_mask]
    print(config['pathways']['type'], config['pathways']['value'], 'X', X_df.shape, 'Y', Y_df.shape, file=log)

if config['DE_filter'] is not None:
    assert re.match('^LFC|^pval|^max_rank', config['DE_filter']['type'])
    if config['DE_filter']['type'].startswith('LFC'):
        top_n = np.argsort(peak_annot_df[config['DE_filter']['type']].abs().values)[::-1]
    else:
        top_n = np.argsort(peak_annot_df[config['DE_filter']['type']].values)
    X_df = X_df.iloc[:, top_n[:config['DE_filter']['value']]]
    print(config['DE_filter']['type'], config['DE_filter']['value'], 'X', X_df.shape, 'Y', Y_df.shape, file=log)

assert X_df.index.equals(Y_df.index) and (donors is None or X_df.index.equals(donors.index))

nested_FS = 'nested_FS' in config and config['nested_FS']
random_state = config['seed']
scale = False #  config['scale']
n_pcs = config['pca']
pca_first = config['pca_first'] if 'pca_first' in config else False
test_splits = config['test_splits']
cv_splits = config['cv_splits']
scoring = config['scoring']
alt_scoring = config['alt_scoring']

if config['select'] is not None:
    if not nested_FS:
        _fs_score = bcg_utils.FEATURE_SELECTION_SCORES[config['select']['type']]
        if config['select']['type'].startswith('spearman_rho_with_bootstrap'):
            feature_selection_score = lambda X, y, k: _fs_score(X, y, random_state)
        elif config['select']['type'].startswith('f_classif_with_low_variance'):
            feature_selection_score = lambda X, y, k: _fs_score(X, y, config['select']['value'], random_state)
        elif config['select']['type'].endswith('_stability'):
            feature_selection_score = lambda X, y, k: _fs_score(X, y, config['select']['value'], random_state)
        elif 'blind_DE' in config['select']['type']:
            assert Y_df.shape[1] == 1
            feature_selection_score = lambda X, y, k: _fs_score(k, random_state if test_splits != -1 else 100,
                                                                Y_df.columns[0], X_df.columns.values)
        else:
            feature_selection_score = lambda X, y, k: _fs_score(X, y)
    else:
        # Cannot be lambda!
        feature_selection_score = bcg_utils.PICKABLE_SCORES[config['select']['type']]

    n_features = config['select']['value']
else:
    feature_selection_score = None
    n_features = None

# C = config['C']
# l1_ratio = config['l1_ratio']
# alpha = config['alpha']
# n_estimators = config['n_estimators']
# learning_rate = config['learning_rate']

_Estimator, est_kwargs, param_names = ESTIMATORS[config['estimator']]
if 'random_state' in est_kwargs:
    est_kwargs['random_state'] = random_state
estimator = _Estimator(**est_kwargs) # max_iter=1000,
param_grid = {p: config[p] for p in param_names}

if isinstance(estimator, RidgeClassifierCV):
    estimator.set_params(alphas=config['alpha'], scoring='accuracy')
    assert cv_splits == -1

if isinstance(estimator, BaggingClassifier):
    estimator.base_estimator.set_params(Cs=config['C'], cv=cv_splits, random_state=random_state)
    print(estimator.base_estimator)

print(estimator)

_KFold = GroupKFold if donors is not None else (StratifiedKFold if bcg_utils.is_classifier(estimator) else SortOfStratifiedKFold if config['scoring'].startswith('quantile_') else KFold)
_LOO = LeaveOneGroupOut if donors is not None else LeaveOneOut
cv_kwargs = dict() if donors is not None else dict(shuffle=True, random_state=None)
# for debugging
# X_df = X_df.iloc[:, np.random.RandomState(random_state).randint(0, high=X_df.shape[1], size=min(X_df.shape[1], 100))]
# print('X_df', X_df.shape, file=log)

# print('Random features')
# X_df = X_df.iloc[:, np.random.RandomState(random_state).randint(0, high=X_df.shape[1], size=min(X_df.shape[1], 189))]

for _the_target in _real_targets if _real_targets is not None else Y_df.columns:
    print('\n----------------\n{}\n----------------\n'.format(_the_target), file=log)
    if _the_target.startswith('IMPUTED_') and _the_target[len('IMPUTED_'):] in Y_df:
        target_runs = [('Check non-imputed samples', _the_target[len('IMPUTED_'):]), ('main', _the_target)]
    else:
        target_runs = [('main', _the_target)]
    assert target_runs[-1][0] == 'main'
    non_imputed_samples, non_imputed_quantiles = None, None
    for run_type, target in target_runs:
        if run_type != 'main':
            print('{}: {}'.format(run_type, target), file=log)
        # if target in ['DONOR:IC_DATE_REAL', 'DONOR:IC_MONTH_REAL', 'DONOR:IC_TIME_REAL', 'SAMPLE:VISIT_DATE_REAL', 'SAMPLE:VISIT_MONTH_REAL', 'SAMPLE:VISIT_TIME_REAL']:
        if config['remove_evening']:
            evening_cohort_removed = True
            if isinstance(Y_df.index, pd.MultiIndex):
                remove_evening = Y_df.index.get_level_values('SAMPLE:DONOR').isin(evening_donors)
            else:
                remove_evening = Y_df.index.isin(evening_donors)
        else:
            evening_cohort_removed = False
            remove_evening = np.zeros((Y_df.shape[0]), dtype=bool)

        if config['remove_300BCG315']:
            BCG315_removed = True
            if isinstance(Y_df.index, pd.MultiIndex):
                BCG315 = Y_df.index.get_level_values('SAMPLE:DONOR').isin(['300BCG315'])
            else:
                BCG315 = Y_df.index.isin(['300BCG315'])
        else:
            BCG315_removed = False
            BCG315 = np.zeros((Y_df.shape[0]), dtype=bool)

        models_fn = fn.format(data='{}.models'.format(target), ext='pickle.gz')
        predictions_fn = fn.format(data='{}.predictions'.format(target), ext='npz')
        if True or not (os.path.exists(models_fn) and os.path.exists(predictions_fn)):
            try:
                if Y_df[target].dtype == float or Y_df[target].dtype == int:
                    assert X_df.index.is_unique and Y_df.index.is_unique
                    assert X_df.index.equals(Y_df.index) and (donors is None or X_df.index.equals(donors.index))

                    keep = (~Y_df[target].isnull().values) & (~X_df.isnull().values.any(axis=1)) & (~remove_evening) & (~BCG315)
                    y_df = Y_df.loc[keep, target].copy()

                    # spearman_rhos = feature_correlation(X=X_df.loc[keep].values, y=Y_df.loc[keep, target].values,
                    #                                     method='spearman', absolute=True, n_bootstrap=0)
                    # ax = sns.distplot(spearman_rhos, kde=False, bins=50)
                    # if n_features is not None:
                    #     ax.axvline(np.sort(spearman_rhos)[-n_features] if n_features < X_df.shape[1] else 0)
                    # ax.set_xlabel('Spearman correlation of {} features\nwith {}'.format(config['X_name'], target))
                    # ax.set_ylabel('Frequency')
                    # sns.despine()
                    # savefig(fn.format(data='{}.target_feature_correlation.distplot'.format(target), ext='pdf'))
                    # if show:
                    #     plt.show()
                    # else:
                    #     plt.close()

                    if (bcg_utils.is_classifier(estimator) and (config['binarize'] or np.unique(y_df).tolist() == [0, 1])) \
                            or (not bcg_utils.is_classifier(estimator) and np.unique(y_df).tolist() not in [[0, 1], [-1, 1]]):
                        if config['binarize']:
                            if config['binarize'].get('thresholds'):
                                assert not y_df.isnull().any()
                                pcts = y_df.rank(pct=True)
                                pct_negative = pcts.loc[y_df < config['binarize']['negative']].max()
                                pct_positive = pcts.loc[y_df >= config['binarize']['positive']].min()
                                thr_negative = y_df.loc[pcts == pct_negative].max()
                                thr_positive = y_df.loc[pcts == pct_positive].min()
                                quantiles = pd.Series([y_df.min(), thr_negative, thr_positive, y_df.max()],
                                                      index=[0, pct_negative, pct_positive, 1], name='binarizing_quantiles')
                                config['binarize']['negative'] = [0, pct_negative]
                                config['binarize']['positive'] = [pct_positive, 1]
                            else:
                                quantiles = bcg_utils.get_binarize_thresholds(
                                    y_df, negative=config['binarize']['negative'],
                                    positive=config['binarize']['positive']).rename('binarizing_quantiles')
                            print(quantiles)
                            assert isinstance(quantiles, pd.Series)

                            if not config['scoring'].startswith('quantile_'):
                                if run_type == 'main' and non_imputed_quantiles is not None:
                                    assert not config['binarize'].get('thresholds')
                                    print('\n{}'.format(quantiles), file=log)
                                    print(non_imputed_quantiles, file=log)
                                    for q in [config['binarize']['negative'][0], config['binarize']['positive'][0]]:
                                        quantiles.loc[q] = min(quantiles.loc[q], non_imputed_quantiles.loc[q])
                                    for q in [config['binarize']['negative'][1], config['binarize']['positive'][1]]:
                                        quantiles.loc[q] = max(quantiles.loc[q], non_imputed_quantiles.loc[q])
                                    print('{}\n'.format(quantiles), file=log)
                                else:
                                    non_imputed_quantiles = quantiles.copy()
                                y_df = bcg_utils.binarize_with_thresholds(y_df, quantiles,
                                                                negative=config['binarize']['negative'],
                                                                positive=config['binarize']['positive'])

                        binarized = ~y_df.isnull().values
                        assert X_df.loc[keep].loc[binarized].index.equals(y_df.loc[binarized].index)

                        # ax = sns.distplot(Y_df.loc[keep, target], kde=False, label='removed', bins=50)
                        # if not binarized.all():
                        #     ax = sns.distplot(Y_df.loc[keep, target].loc[binarized], kde=False, label='binarized', bins=50)
                        #     ax.legend()
                        # ax.set_xlabel(target)
                        # ax.set_ylabel('Frequency')
                        # sns.despine()
                        # bcg_utils.savefig(fn.format(data='{}.distplot'.format(target), ext='pdf'))
                        # if show:
                        #     plt.show()
                        # else:
                        #     plt.close()

                        X = X_df.loc[keep].loc[binarized].copy()
                        y = y_df.loc[binarized].copy()
                        assert X.index.equals(y.index)
                        y_real = Y_df.loc[keep, target].loc[binarized].copy()
                        assert y.index.equals(y_real.index)
                        X_extra = X_df.loc[keep].loc[~binarized].copy()
                        y_extra = Y_df.loc[keep, target].loc[~binarized].copy()
                        assert X_extra.index.equals(y_extra.index)
                        assert not X_extra.index.isin(X.index).any()

                        groups = donors.loc[keep].loc[binarized].copy() if donors is not None else None
                        assert X.index.equals(y.index) \
                               and (groups is None or (groups.index.equals(y.index) and not groups.isnull().any()))
                        if run_type == 'main':
                            if evening_cohort_removed:
                                print('Removed evening cohort: {} samples'.format(remove_evening.sum()), file=log)
                            if BCG315_removed:
                                print('Removed 300BCG315: {} samples'.format(BCG315.sum()), file=log)
                            print(target, 'X', X.shape, 'y', y.shape, 'groups', groups.shape if groups is not None else None)
                            print(target, 'X', X.shape, 'y', y.shape, 'groups', groups.shape if groups is not None else None, file=log)

                        if X.shape[0] > 1 and np.unique(y).shape[0] > 1:
                            if run_type == 'main':
                                print('{}{} {}'.format(estimator.__class__.__name__,
                                                   '_{}'.format(estimator.penalty) if hasattr(estimator, 'penalty') else '',
                                                   param_grid), file=log)

                            shuffle = np.random.RandomState(random_state if test_splits != -1 else 100).permutation(X.shape[0])
                            X, y, groups = X.iloc[shuffle], y.iloc[shuffle], groups.iloc[shuffle] if groups is not None else None

                            if bcg_utils.is_classifier(estimator):
                                assert all((y == 1) | (y == 0))
                                majority_class = bcg_utils.get_majority_class(y.values)
                                majority_class_mask = (y == majority_class).values
                                majority_class_count = majority_class_mask.sum()
                                assert majority_class_count >= y.shape[0] / 2
                                if run_type == 'main':
                                    print('Majority class {:.0f}% ({} vs. {})'.format(round((majority_class_count / y.shape[0]) * 100), majority_class_count, y.shape[0] - majority_class_count))
                                    print('Majority class {:.0f}% ({} vs. {})'.format(round((majority_class_count / y.shape[0]) * 100), majority_class_count, y.shape[0] - majority_class_count), file=log)
                                if config['downsampling'] and majority_class_count > y.shape[0] / 2:
                                    if run_type == 'main':
                                        print('Before downsampling {} ({}/{})'.format(
                                            y.shape[0], majority_class_count, y.shape[0] - majority_class_count), file=log)
                                    remove_n = majority_class_count - (y.shape[0] - majority_class_count)
                                    X = pd.concat([X.loc[majority_class_mask].iloc[remove_n:], X.loc[~majority_class_mask]])
                                    y = pd.concat([y.loc[majority_class_mask].iloc[remove_n:], y.loc[~majority_class_mask]])
                                    if groups is not None:
                                        groups = pd.concat([groups.loc[majority_class_mask].iloc[remove_n:], groups.loc[~majority_class_mask]])
                                    majority_class_count = (y == majority_class).sum()
                                    if run_type == 'main':
                                        print('After downsampling {} ({}/{})'.format(y.shape[0], majority_class_count,
                                                                                 y.shape[0] - majority_class_count), file=log)

                            kfold_kwargs = dict(thresholds=quantiles.drop([0, 1]).values) if config['scoring'].startswith('quantile_') else dict()
                            n_smallest_class = np.min([(y == c).sum() for c in np.unique(y)]) if bcg_utils.is_classifier(estimator) else len(y)
                            test_kfold = _KFold(n_splits=min(test_splits, n_smallest_class), **kfold_kwargs) if test_splits != -1 else _LOO()
                            if not bcg_utils.is_classifier(estimator) and len(y) == test_splits + 1:
                                _corrected_splits = len(y) - 2
                            else:
                                _corrected_splits = n_smallest_class - 1
                            cv_kfold = _KFold(n_splits=min(cv_splits, _corrected_splits), **cv_kwargs, **kfold_kwargs) if cv_splits != -1 else _LOO()
                            if run_type == 'main':
                                print('test_kfold:', test_kfold, file=log)
                                print('cv_kfold:', cv_kfold, file=log)
                            # approx_train_fold_size = np.floor(np.floor(X.shape[0] / test_splits * (test_splits - 1)) / cv_splits * (cv_splits - 1))

                            assert X.index.equals(y.index) and (groups is None or (groups.index.equals(y.index) and not groups.isnull().any()))

                            if config.get('bootstrap'):
                                print('Bootstrap:', X.shape, 'unique:', len(X.index.str.split(':', expand=True).get_level_values(0).drop_duplicates()) / len(X))

                            # if config.get('subsample'):
                            #     if bcg_utils.is_classifier(estimator):
                            #         assert set(y) in [set([0, 1]), set([-1, 1])]
                            #         _idx1 = X.loc[y == 1].sample(frac=config.get('subsample'), replace=False, random_state=config['seed']).index
                            #         _idx2 = X.loc[y != 1].sample(frac=config.get('subsample'), replace=False, random_state=config['seed']).index
                            #         assert len(_idx1.intersection(_idx2)) == 0
                            #         _idx = np.concatenate([_idx1, _idx2])
                            #     else:
                            #         _idx = X.sample(frac=config.get('subsample'), replace=False, random_state=config['seed']).index
                            #     print('Subsampling', X.shape, '-->', end=' ')
                            #     X, y, y_real, groups = X.loc[_idx], y.loc[_idx], y_real.loc[_idx], groups.loc[_idx] if groups is not None else None
                            #     print(X.shape)

                            assert X.index.equals(y.index) and (groups is None or (groups.index.equals(y.index) and not groups.isnull().any()))

                            if run_type == 'main':
                                log.flush()
                                grids, y_true, y_pred, samples_pred, feature_importances, selected_features, train_scores, test_scores = \
                                    _nested_cv(nested_FS=nested_FS,
                                               estimator=estimator,
                                               X=X.values,
                                               y=y.values,
                                               y_real=y_real.values,
                                               X_extra=X_extra.values,
                                               y_extra=y_extra.values,
                                               groups=groups.values if groups is not None else None,
                                               samples=X.index.values,
                                               peaks=X.columns.values,
                                               target=target,
                                               test_kfold=test_kfold, cv_kfold=cv_kfold,
                                               random_state=random_state,
                                               param_grid=param_grid,
                                               scoring=scoring,
                                               alt_scoring=alt_scoring,
                                               permute_labels=config.get('permute_labels'),
                                               subsample=config.get('subsample'),
                                               scale=scale,
                                               binarized=quantiles.drop([0, 1]) if scoring.startswith('quantile_') else None,
                                               merged_folds_scoring=config['merge_folds'] if 'merge_folds' in config else False,
                                               feature_selection_name=config['select']['type'] if config['select'] is not None else None,
                                               feature_selection_mask=X.columns.str.startswith('CONS0'),
                                               feature_selection_score=feature_selection_score,
                                               n_features=n_features,
                                               random_features=config['select']['type'].startswith('rand_DE_') if config['select'] is not None else False,
                                               # n_pcs=np.min(approx_train_fold_size, X.shape[1], n_pcs),
                                               n_pcs=n_pcs, pca_first=pca_first,
                                               n_jobs=n_jobs, verbose=1, out=log)

                                feature_names = []
                                for i in range(len(grids)):
                                    if n_pcs is None:
                                        if n_features is None:
                                            feature_names.append(X_df.columns.values)
                                        else:
                                            feature_names.append(X_df.columns.values[selected_features[i]])
                                    else:
                                        feature_names.append(['PC{}'.format(i + 1) for i in range(
                                            grids[i].best_estimator_.named_steps.pca.n_components_
                                        )])

                                np.savez(predictions_fn,
                                         y_true=y_true,
                                         y_pred=y_pred,
                                         samples=samples_pred,
                                         non_imputed_samples=non_imputed_samples,
                                         feature_importances=feature_importances,
                                         feature_names=feature_names,
                                         train_scores=train_scores,
                                         test_scores=test_scores)
                                with gzip.open(models_fn, 'wb') as f:
                                    joblib.dump(grids, f)
                            else:
                                non_imputed_samples = X.index.values
                                print('There are {} non-imputed samples'.format(len(non_imputed_samples)), file=log)
                        else:
                            if run_type == 'main':
                                print('Skipping {} (number of samples: {}, number of unique labels: {})'.format(target, X.shape[0], np.unique(y).shape[0]), file=log)
                            else:
                                non_imputed_samples = np.asarray([])
                                print('There are {} non-imputed samples'.format(len(non_imputed_samples)), file=log)
                    else:
                        if run_type == 'main':
                            print('Skipping {} because it is a binary target'.format(target), file=log)
            except Exception as e:
                print('\nERROR', file=log)
                traceback.print_exception(*sys.exc_info(), file=log)
                raise e
        else:
            if run_type == 'main':
                print(models_fn)
                print(predictions_fn)
                print('Skipping {} because it has already been processed'.format(target), file=log)

print('\nDone.', file=log)
log.flush()
log.close()
print('Done.')
