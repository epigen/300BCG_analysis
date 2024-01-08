import yaml
import bcg_utils
from misc import *

FEATURE_SELECTION = [dict(type=t, value=v) for t in [
    'random'
] for v in [100, 500, 1000, 2000, 5000]]
NESTED_FS = True

name = '{model}' + os.path.sep + '{celltype}' + os.path.sep + '{X}_{X_visits}:{Y}_{Y_visits}'
C = [1e6, 1e5, 1e4, 1e3, 1e2, 1e1, 1e0, 1e-1, 1e-2, 1e-3, 1e-4, 1e-5, 1e-6]
# l1_ratio = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9]
l1_ratio = [0.5]

degree = [2]
n_estimators = [1000]
learning_rate = [0.01, 0.05, 0.1, 0.2]
epsilon = [1.05, 1.1, 1.15, 1.2, 1.25, 1.3, 1.35, 1.4, 1.45, 1.5]

# UMAP default is 15 and 0.1
UMAP_n_neighbors = 15
UMAP_min_dist = 0.1

DONOR_REGEX = r'^DONOR:IC_TIME_REAL$|^DONOR:AGE$|^DONOR:BMI$|^DONOR:SEX$|^DONOR:oralContraceptivesIncludingMen$|^DONOR:IC_alcoholInLast24h$'
SEASON_REGEX = r'^DONOR:IC_DATE_2PI_SIN$|^DONOR:IC_DATE_2PI_COS$'
q1, q2 = 0.25, 0.75
celltype = 'PBMC'
bootstrap = False
merge_folds = True
test_splits = -1
cv_splits = -1
group_by_date = -1

for model, permute_labels, subsample, n_seeds in [
    ('final_V1_corrected_combat_groupByExactDate_thr', False, False, 1),
    ('final_V1_corrected_combat_groupByExactDate_thr', False, 0.9, 100),
    ('final_V1_corrected_combat_groupByExactDate_thr', True, False, 100),
    ('final_V1_corrected_combat_groupByExactDate_thr', True, 0.9, 100),
]:
    model += '_subsample{}'.format(subsample) if subsample else ''
    model += '_permuteLabels' if permute_labels else ''
    for est in ['LR_L2']:  # 'LR_L2', 'RFC', 'LR_Enet', 'SVM'
        if est == 'LR_Enet' and l1_ratio == [0.5]:
            C = [4 * c for c in [1e6, 1e5, 1e4, 1e3, 1e2, 1e1, 1e0, 1e-1, 1e-2, 1e-3, 1e-4, 1e-5, 1e-6]]
        for X_visits, Y_visits in [
            (['V1'], ['V1']),
            # (['V3-V1'], ['V1']),
        ]:
            config_list_fn = os.path.join(ML_RESULTS_ROOT, 'ML_{}_{}_{}_{}'.format(
                model, ''.join(X_visits), ''.join(Y_visits), est))
            print(config_list_fn)
            with open(config_list_fn, 'w') as todo_file:
                i = 0
                for use_thresholds, q1, q2 in [(True, 1.2, 1.6), (True, 1.2, 1.7), (True, 1.2, 1.8)]:
                    if use_thresholds:
                        quartiles = dict(negative=float(np.log2(q1)), positive=float(np.log2(q2)), thresholds=True)
                    else:
                        quartiles = dict(negative=[0.0, q1], positive=[q2, 1.0])
                    for Y in ['thm.innate_nonspecific_24h_wo_LAC_IL10_IL1ra_V3'] + \
                             (['thm.adaptive_MTB_7d_V3'] if not subsample else []):
                        _classification = Y in ['DONOR:SEX', 'SAMPLE:VISIT'] \
                                          or Y.startswith('CYTO:') \
                                          or Y.startswith('IMPUTED_CYTO:') \
                                          or Y.startswith('CM:') \
                                          or Y.startswith('thm.')
                        assert _classification
                        for X_name, X, scale, dim_reduce, n_components, poly in [
                            ('ATAC_SEASON_DONOR', ['ATAC', SEASON_REGEX, DONOR_REGEX], True, 'PCA', 5, None),
                        ] + (
                                [
                                    ('DONOR', [DONOR_REGEX], True, None, None, None),
                                    ('DONOR_SEASON', [DONOR_REGEX, SEASON_REGEX], True, None, None, None),
                                ]
                                if not permute_labels else
                                []
                        ):
                            if Y.startswith(X_name) and X_visits == Y_visits:
                                continue

                            for seed in range(n_seeds):
                                seed = (seed + 1) * 100
                                promoters = []
                                _pathways = [None]  # [dict(type='GO_glycolysis_gluconeogenesis_regions', value='GO_glycolysis_gluconeogenesis_regions')] # [None]
                                DE_regions = [None] if X_name == 'ATAC' else [None]
                                feature_selection = [None]  # FEATURE_SELECTION if X_name == 'ATAC' else [None]
                                for features in [None] + promoters:
                                    for pathways in _pathways:
                                        for DE_filter in DE_regions:
                                            for select in [None] if DE_filter else feature_selection:  # dict(type='f_classif_stability', value=10)  # f_classif_stability
                                                if _classification:
                                                    C, scoring, alt_scoring, downsampling = \
                                                        C, 'roc_auc', 'average_precision', False  # partial_roc_auc
                                                    estimators = [est]
                                                    binarize = None if Y in get_basic_binary_specials().keys() else quartiles
                                                else:
                                                    C, scoring, alt_scoring, downsampling = \
                                                        C, 'quantile_roc_auc', 'quantile_average_precision', False
                                                    assert (scoring.startswith('quantile_') and alt_scoring.startswith('quantile_')) or (not scoring.startswith('quantile_') and not alt_scoring.startswith('quantile_'))
                                                    estimators = [est]
                                                    binarize = quartiles if scoring.startswith('quantile_') else None

                                                for estimator in estimators:
                                                    config = {}
                                                    config['name'] = name
                                                    config['celltype'] = celltype
                                                    config['model'] = model
                                                    config['X_name'] = X_name
                                                    config['X'] = X
                                                    config['Y'] = '^{}$'.format(Y)
                                                    config['X_visits'] = X_visits
                                                    config['Y_visits'] = Y_visits
                                                    config['bootstrap'] = bootstrap
                                                    config['downsampling'] = downsampling
                                                    config['binarize'] = binarize
                                                    config['features'] = features
                                                    config['pathways'] = pathways
                                                    config['DE_filter'] = DE_filter
                                                    config['seed'] = seed
                                                    config['scale'] = scale
                                                    config['pca_first'] = True
                                                    config['group_by_date'] = group_by_date
                                                    config['merge_folds'] = merge_folds
                                                    config['poly'] = poly
                                                    config['pca'] = n_components if dim_reduce is None else None
                                                    config['global_KPCA_degree'] = 2
                                                    config['global_KPCA'] = n_components if dim_reduce == 'KPCA' else None
                                                    config['global_PCA'] = n_components if dim_reduce == 'PCA' else None
                                                    config['global_UMAP'] = n_components if dim_reduce == 'UMAP' else None
                                                    config['n_neighbors'] = UMAP_n_neighbors if dim_reduce == 'UMAP' else None
                                                    config['min_dist'] = UMAP_min_dist if dim_reduce == 'UMAP' else None
                                                    config['select'] = select
                                                    config['test_splits'] = test_splits
                                                    config['cv_splits'] = cv_splits
                                                    config['scoring'] = scoring
                                                    config['alt_scoring'] = alt_scoring
                                                    config['estimator'] = estimator
                                                    config['C'] = sorted(C)
                                                    config['alpha'] = sorted(C)[::-1]
                                                    config['l1_ratio'] = l1_ratio
                                                    config['degree'] = degree
                                                    config['n_estimators'] = n_estimators
                                                    config['learning_rate'] = learning_rate
                                                    config['epsilon'] = epsilon
                                                    config['subsample'] = subsample
                                                    config['permute_labels'] = permute_labels
                                                    config['remove_evening'] = (Y == 'SAMPLE:IC_TIME_REAL')
                                                    config['remove_300BCG315'] = False
                                                    config['nested_FS'] = NESTED_FS

                                                    run_dir = bcg_utils.make_dir(ML_RESULTS_ROOT, *get_data_name_from_config(config).split(os.path.sep))
                                                    fn = get_filename_from_config(config, run_dir).format(data='config', ext='yml')
                                                    with open(fn, 'w') as f:
                                                        yaml.dump(config, f, default_flow_style=False, sort_keys=False)
                                                    print(fn, file=todo_file)
                                                    i += 1
        print(i)
