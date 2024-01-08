from __future__ import print_function
from __future__ import absolute_import
from __future__ import division

from abc import ABC

from six import itervalues, iteritems, string_types
from builtins import int, str, bytes, chr, object, dict
from builtins import range, map, zip, filter

import traceback
import re
from sklearn.linear_model import RidgeClassifierCV
from sklearn.ensemble import BaggingClassifier
import sys
import yaml
import numpy as np
from scipy import stats
from statsmodels.stats.multitest import multipletests
import pandas as pd
from pandas.api.types import is_bool_dtype, is_integer_dtype, is_string_dtype, is_float_dtype
import itertools
import os
import gzip
import pickle
import collections
import copy
import time
import seaborn as sns
import sklearn
from functools import wraps
from statsmodels.formula.api import ols, mixedlm
from statsmodels.stats.anova import anova_lm
from collections import OrderedDict
from collections import Sized, defaultdict
from functools import partial
from math import sqrt
import matplotlib.gridspec
from concurrent.futures import ProcessPoolExecutor, as_completed
from contextlib import contextmanager
from collections import Iterable
import typing as t
from itertools import takewhile

import scipy.sparse as sp

from sklearn.base import BaseEstimator
# from sklearn.feature_selection import SelectorMixin
from sklearn.utils import check_array
from sklearn.utils.validation import check_is_fitted
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import label_binarize
from sklearn.metrics import roc_auc_score, r2_score, accuracy_score, mean_squared_error, average_precision_score, auc
from sklearn.feature_selection import VarianceThreshold
from sklearn.model_selection import KFold, StratifiedKFold, GroupKFold, LeaveOneGroupOut, LeaveOneOut
from sklearn.feature_selection import mutual_info_classif, mutual_info_regression, f_classif, f_regression
from sklearn.model_selection import GridSearchCV
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
from sklearn.pipeline import Pipeline
from sklearn.metrics import get_scorer
from sklearn.metrics import make_scorer
import joblib
from statsmodels.tools.sm_exceptions import ConvergenceWarning
import warnings
from sklearn.feature_selection import SelectKBest
from sklearn.model_selection import cross_val_predict as _cross_val_predict
from sklearn.model_selection._validation import _check_is_permutation
from sklearn.model_selection._validation import _num_samples
from sklearn.preprocessing import LabelEncoder
from sklearn.model_selection import learning_curve
from sklearn.model_selection import GridSearchCV
# from sklearn.model_selection._search import _check_param_grid
from sklearn.model_selection._split import check_cv
from sklearn.base import is_classifier, clone
from sklearn.utils.validation import indexable
from sklearn.model_selection._validation import _fit_and_score, _fit_and_predict
#from sklearn.utils.fixes import rankdata
#from sklearn.utils.fixes import MaskedArray
from sklearn.metrics import check_scoring
from sklearn.model_selection import ParameterGrid
from sklearn.preprocessing import StandardScaler
from sklearn.linear_model import LogisticRegression
from sklearn.linear_model import ElasticNet
from sklearn.linear_model import Lasso
from sklearn.linear_model import lasso_path
from sklearn.linear_model import SGDClassifier
from sklearn.linear_model import Ridge, RidgeClassifier, LinearRegression
from sklearn.kernel_ridge import KernelRidge
from sklearn.gaussian_process import GaussianProcessClassifier, GaussianProcessRegressor
from sklearn.naive_bayes import GaussianNB, MultinomialNB, BernoulliNB
from sklearn.svm import SVC
from sklearn.metrics.pairwise import euclidean_distances
from sklearn.metrics._scorer import _BaseScorer
from sklearn.ensemble import RandomForestRegressor, RandomForestClassifier
from sklearn.neighbors import KNeighborsRegressor, KNeighborsClassifier
from sklearn.decomposition import PCA, KernelPCA
from sklearn.base import BaseEstimator, TransformerMixin
from sklearn.preprocessing import LabelBinarizer, MultiLabelBinarizer
from matplotlib.colors import ListedColormap
import matplotlib.pyplot as plt
import matplotlib as mpl
from tempfile import TemporaryDirectory
import subprocess
from bcg_colors import BCG_PAL, BLUE, RED, LIGHT_GRAY, PROMOTER_MAP_COL, DISTAL_MAP_COL
import warnings

####################
# scipy.stats
####################
import scipy
from scipy.stats import boxcox, rankdata
from scipy.stats import pearsonr, spearmanr
from scipy.stats import shapiro
from scipy.stats import ttest_ind, mannwhitneyu, f_oneway, kruskal, ttest_rel, wilcoxon, chi2_contingency, fisher_exact

# pearsonr(a, b)
# spearmanr(a, b, nan_policy='omit')

# Independent samples (two groups):
# scipy.stats.ttest_ind(a, b, nan_policy='omit') # parametric
# scipy.stats.mannwhitneyu(a, b) # non-parametric

# Independent samples (more than two groups):
# scipy.stats.f_oneway(x, y, z) # parametric
# scipy.stats.kruskal(x, y, z, nan_policy='omit') # non-parametric

# Related samples:
# scipy.stats.ttest_rel(a, b, nan_policy='omit') # parametric
# scipy.stats.wilcoxon(a, b) # non-parametric


def _nanrankdata(a, axis=-1, inplace=False):

    if hasattr(a, "dtype") and np.issubdtype(a.dtype, np.integer):
        raise ValueError("Integer type is not supported.")

    if isinstance(a, (tuple, list)):
        if inplace:
            raise ValueError("Can't use `inplace=True` for {}.".format(type(a)))
        a = np.asarray(a, float)

    orig_shape = a.shape
    if a.ndim == 1:
        a = a.reshape(orig_shape + (1,))

    if not inplace:
        a = a.copy()

    def rank1d(x):
        idx = ~np.isnan(x)
        x[idx] = rankdata(x[idx])
        return x

    a = a.swapaxes(1, axis)
    a = np.apply_along_axis(rank1d, 0, a)
    a = a.swapaxes(1, axis)

    return a.reshape(orig_shape)


def quantile_gaussianize(X, axis=1):

    X = X.copy()
    orig_shape = X.shape
    if X.ndim == 1:
        X = X.reshape(orig_shape + (1,))

    D = X.swapaxes(1, axis)
    D = np.ma.masked_invalid(D)
    D *= -1

    with warnings.catch_warnings():
        warnings.filterwarnings("ignore", category=RuntimeWarning)
        D = _nanrankdata(D)
        D = D / (np.isfinite(D).sum(axis=0) + 1)
        D = np.apply_along_axis(scipy.stats.norm.isf, 0, D)

    D = D.swapaxes(1, axis)
    X[:] = D

    return X.reshape(orig_shape)


def omit_nans(f, dependent, verbose=True):
    """
    pearsonr = omit_nans(pearsonr, dependent=True)
    wilcoxon = omit_nans(wilcoxon, dependent=True)
    f_oneway = omit_nans(f_oneway, dependent=False)
    """
    @wraps(f)
    def wrapper(*args):
        data = [np.asarray(d) for d in takewhile(lambda x: isinstance(x, Iterable), args)]
        print('{} input arrays'.format(len(data)))
        if np.isnan(data).any():
            if verbose:
                print('Removing nans from {} input arrays'.format(len(data)))
            assert len(data) != 0
            assert np.all([len(d.shape) == 1 for d in data])
            if not dependent:
                data = [d[~np.isnan(d)] for d in data]
            else:
                assert len(data) > 1
                assert np.all([data[0].shape[0] == d.shape[0] for d in data])
                not_nan = np.ones_like(data[0], dtype=bool)
                for d in data:
                    not_nan &= ~np.isnan(d)
                data = [d[not_nan] for d in data]
            return f(*data, *args[len(data):])
        else:
            return f(*args)

    return wrapper

####################

CONCURRENT = False

OLS_REPLACE_CHARS = ':./'

BLUE_GREEN_CMAP = sns.cubehelix_palette(start=.5, rot=-.75, light=.75, as_cmap=True)
PLOT_CONTEXTS = dict(paper=0.8, notebook=1, talk=1.5, poster=2)

ENRICHR_COLS = ['Gene_set', 'Term', 'P-value', 'Adjusted P-value', 'Combined Score', 'Odds Ratio', 'Overlap', 'Genes']



PALETTE_CODES = {
    'diverging': 'RdBu_r',
    'sequential': 'cubehelix',
    'qualitative': BCG_PAL,

    'D': 'RdBu_r',
    'S': 'cubehelix',
    'Q': BCG_PAL
}

R, NR = 'Responder', 'Non-responder'
CR, PR, SD, PD = 'CR', 'PR', 'SD', 'PD'
RECIST_RESPONSES = [CR, PR, SD, PD]
BINARY_RESPONSES = [R, NR]
BORTEZOMIB_RESPONSES = [CR, PR, 'MR', 'NC', PD]

EN_DASH = u'\u2013'
EM_DASH = u'\u2014'
LESS_EQUAL = u'\u2264'
GREATER_EQUAL = u'\u2265'
MUCH_GREATER = u'\u226B'
MUCH_LESS = u'\u226A'
NOT_EQUAL = u'\u2260'
ARROW = u'\u2192'
UP_ARROW = u'\u2191'
DOWN_ARROW = u'\u2193'
DOUBLE_ARROW = u'\u2194'
DELTA = u'\u0394'
SMALL_DELTA = u'\u03B4'
PLUS_MINUS = u'\u00B1'
UNION = u'\u222A'
MINUS = u'\u2212'
PLUS = u'\u002B'
ALPHA = u'\u03B1'
BETA = u'\u03B2'
GAMMA = u'\u03B3'
NBSPC = u'\u00A0'
PI = u'\u03C0'
TIMES = u'\u00D7'

EST_STEP = 'estimator'
ESTIMATOR_PREFIX = '{}__'.format(EST_STEP)
PCA_STEP = 'pca'
SCALE_STEP = 'scaler'
FS_STEP = 'feature_selector'

TSS_PROXIMAL = 'TSS_PROXIMAL'
PROMOTERS = 'PROMOTERS'
GENE_AND_DISTAL_10kb = 'GENE_AND_DISTAL_10kb'
DISTAL_1Mb = 'DISTAL_1Mb'
ALL_GENES = 'single_gene_1Mb_wo_manual'

TSS_PROXIMAL_GROUP = ['TSS', 'TSS_proximal', 'TSS_overlap', 'TSS_FIP']
GENE_BODY_GROUP = ['gene_body', 'gene_body_FIP', 'gene_body_overlap_end', 'gene_body_overlap_start']
DISTAL_GROUP = ['distal']
INTERGENIC_GROUP = ['intergenic']

PEAKS_TO_GENES = {
    TSS_PROXIMAL: dict(
        location=TSS_PROXIMAL_GROUP,
        or_distance=None,
        and_feat_type=['gene:protein_coding', 'transcript:protein_coding']
    ),
    # PROMOTERS: dict(location=TSS_PROXIMAL_GROUP,
    #                    or_distance=None,
    #                    and_feat_type=['gene:protein_coding', 'transcript:protein_coding'],
    #                    and_reg_feature=['promoter']),
    GENE_AND_DISTAL_10kb: dict(
        location=TSS_PROXIMAL_GROUP +  GENE_BODY_GROUP + DISTAL_GROUP,
        or_distance=None,
        and_feat_type=['gene:protein_coding', 'transcript:protein_coding']
    ),
    # ALL_GENES: dict(location=None,
    #                 or_distance=None,
    #                 and_feat_type=['gene:protein_coding', 'transcript:protein_coding']),
}

LNC_RNA_TSS_PROXIMAL = 'LNC_RNA_TSS_PROXIMAL'
LNC_RNA_GENE_AND_DISTAL_10kb = 'LNC_RNA_GENE_AND_DISTAL_10kb'
LNC_RNA_DISTAL_1Mb = 'LNC_RNA_DISTAL_1Mb'

PEAKS_TO_LNC_RNA = {
    LNC_RNA_TSS_PROXIMAL: dict(
        location=['TSS', 'TSS_proximal'],
        or_distance=None,
        and_feat_type=['gene:lncRNA', 'transcript:lncRNA']),
    LNC_RNA_GENE_AND_DISTAL_10kb: dict(
        location=['TSS', 'TSS_proximal', 'gene_body'],
        or_distance=1e4,
        and_feat_type=['gene:lncRNA', 'transcript:lncRNA']),
    LNC_RNA_DISTAL_1Mb: dict(
        location=['TSS', 'TSS_proximal'],
        or_distance=1e6,
        and_feat_type=['gene:lncRNA', 'transcript:lncRNA']),
}


def peaks_to_genes(df, location=None, or_distance=None, and_feat_type=None, and_reg_feature=None):
    if isinstance(and_feat_type, str):
        and_feat_type = [and_feat_type]
    if isinstance(and_reg_feature, str):
        and_reg_feature = [and_reg_feature]
    if isinstance(location, str):
        location = [location]

    mask = np.ones((df.shape[0],), dtype=bool)
    if location:
        mask &= df['characterization'].isin(location)
    if or_distance:
        mask |= (df['distance'] < or_distance)
    if and_feat_type:
        mask &= df['feat_type'].isin(and_feat_type)
    if and_reg_feature:
        mask &= df['reg_feature'].isin(and_reg_feature)
    return mask


def pearson_r_score(x, y):
    return pearsonr(x, y)[0]


def spearman_r_score(x, y):
    return spearmanr(x, y)[0]


def inverted_spearman_r_score(x, y):
    return spearman_r_score(x, -1 * y)


def make_int_classes(y, n_classes=None):
    if isinstance(y, (list, tuple)):
        y = np.array(y)
    class_labels = np.unique(y)
    assert n_classes is None or len(class_labels) == n_classes
    for i, label in enumerate(class_labels):
        y[y == label] = i
    return y.astype(int)


def auroc_score(y_true, y_score):
    # Convert binary labels of any kind to 0/1 numerical labels
    y_true_num = make_int_classes(y_true, n_classes=2)
    return roc_auc_score(y_true=y_true_num, y_score=y_score)


def abs_auroc_score(y_true, y_score):
    auroc = auroc_score(y_true, y_score)
    return auroc if auroc >= 0.5 else (1 - auroc)


def inverted_auroc_score(y_true, y_score):
    return auroc_score(y_true, -1 * y_score)


def mann_whitney_u_for_classification(y_true, y_pred):
    classes = np.unique(y_true)
    assert len(classes) == 2
    preds_by_class = [None, None]
    for i, label in enumerate(classes):
        preds_by_class[i] = y_pred[y_true == label]
    return mannwhitneyu(preds_by_class[0], preds_by_class[1])[1]


def get_median_trick_gamma(X, Y=None, extended=False):
    if Y is None:
        Y = X
    distances = euclidean_distances(X, Y)
    squared_distances = distances.flatten() ** 2
    if extended:
        gamma = 1.0 / np.percentile(squared_distances, [10, 50, 90])
        assert gamma[1] == 1.0 / np.median(squared_distances), (gamma, 1.0 / np.median(squared_distances))
    else:
        gamma = 1.0 / np.median(squared_distances)
    return gamma


def get_kfold(n_splits, groups=None, stratify=False, adjust_n_splits_for_y=None, shuffle=False, random_state=None):
    assert groups is None or (not stratify and not shuffle and random_state is None)
    if adjust_n_splits_for_y is not None:
        classes = np.unique(adjust_n_splits_for_y)
        n_smaller_class = np.min([np.sum(adjust_n_splits_for_y == label) for label in classes])
        if n_smaller_class < n_splits:
            print('Warning: n_splits ({}) > n_smaller_class ({}), setting n_splits to {}'.format(n_splits, n_smaller_class, n_smaller_class))
            n_splits = n_smaller_class
    _KFfold = GroupKFold if groups is not None else (StratifiedKFold if stratify else KFold)
    _kwargs = dict() if groups is not None else dict(shuffle=shuffle, random_state=random_state)
    return _KFfold(n_splits=n_splits, **_kwargs)


def open_file(filename, mode='r', compresslevel=9):
    assert mode in ['r', 'r+b', 'a', 'a+b', 'w', 'w+b', 'b']
    if filename.endswith('.gz') or (mode in ['r', 'r+b'] and not os.path.exists(filename) and os.path.exists(filename + '.gz')):
        #gzip automatically adds 'b' to the 'r', 'a', and 'w' modes
        return gzip.open(filename if filename.endswith('.gz') else filename + '.gz', mode, compresslevel)
    else:
       return open(filename, mode)


def maybe_gz(fn):
    if not fn.endswith('.gz') and not os.path.exists(fn) and os.path.exists('{}.gz'.format(fn)):
        return '{}.gz'.format(fn)
    else:
        return fn


def line_split_gen(filename, delim='\t', strip='\n', comment='#', skip_rows=0, skip_columns=0):
    with open_file(filename) as f:
        for _ in range(skip_rows):
            f.readline()
        for line in f:
            if comment is not None and line.startswith(comment):
                continue
            line_split = line.strip(strip).split(delim)
            yield line_split[skip_columns:]


def array_upper(list_of_str):
    return array_change_case(list_of_str, case='upper')


def array_lower(list_of_str):
    return array_change_case(list_of_str, case='lower')


def array_change_case(list_of_str, case):
    assert case in ['upper', 'lower']
    return np.array([(s.upper() if case == 'upper' else s.lower()) if not pd.isnull(s) else s for s in list_of_str])


def intersect_index(*lists):
    assert len(lists) > 1
    assert np.all([no_nulls(l) for l in lists])
    assert np.all([is_unique(l) for l in lists])

    merged_df = None
    for i, l in enumerate(lists):
        l_df = pd.DataFrame(np.stack([np.arange(len(l)), l], axis=1), columns=['index{}'.format(i), 'value'])
        merged_df = pd.merge(merged_df, l_df, on='value', how='inner') if merged_df is not None else l_df.copy()

    perm_indexes = [merged_df['index{}'.format(i)].as_matrix().astype(np.int) for i in range(len(lists))]

    assert np.all([np.array_equal(lists[0][perm_indexes[0]], lists[i][perm_indexes[i]]) for i in range(1, len(lists))])

    return perm_indexes


def merge_dicts(old, new):
    d = old.copy()
    d.update(new)
    return d


def binarize_labels_old(y, threshold=None, mean=False, median=False, quartile=False, smaller_class=0):
    assert np.sum([threshold is not None, mean, median, quartile]) == 1
    assert smaller_class in [0, 1]

    y = np.array(y)
    y_binary = np.ones(y.shape) * np.nan
    not_nan = ~np.isnan(y)

    if threshold is not None or mean or median:
        if mean:
            threshold = np.mean(y[not_nan])
        elif median:
            threshold = np.median(y[not_nan])
        y_binary[not_nan] = (y[not_nan] >= threshold if smaller_class == 0 else y[not_nan] < threshold)

    elif quartile:
        q1 = np.percentile(y[not_nan], 25)
        q3 = np.percentile(y[not_nan], 75)
        y_binary_not_nan = y_binary[not_nan]
        y_binary_not_nan[y[not_nan] >= q3] = 1 if smaller_class == 0 else 0
        y_binary_not_nan[y[not_nan] <= q1] = 0 if smaller_class == 0 else 1
        y_binary[not_nan] = y_binary_not_nan

    return y_binary.astype(float)


@contextmanager
def none_context_manager(*args, **kwargs):
    yield None


def _execute(fn, dict_of_args):
    return fn(**dict_of_args)


def maybe_parallel(parallel, n_jobs, fns, dicts_of_args):
    if parallel:
        with ProcessPoolExecutor(max_workers=n_jobs) as executor:
            results = executor.map(_execute, fns, dicts_of_args)
    else:
        results = map(_execute, fns, dicts_of_args)

    return results


def start_timing():
    return time.time()


def elapsed_time(start):
    return time.time() - start


def scp_copy(source, destination):
    os.system('scp {} {}'.format(source, destination))


def load_pickle(filename):
    with open_file(filename, 'r+b') as f:
        return pickle.load(file=f)


def dump_pickle(filename, data, protocol=None, scp=None):
    with open_file(filename, 'w+b') as f:
        pickle.dump(obj=data, file=f, protocol=protocol)
    if scp is not None:
        scp_copy(source=filename, destination=scp)


def savez(filename, scp=None, **kwargs):
    np.savez(file=filename, **kwargs)
    if scp is not None:
        scp_copy(source=filename, destination=scp)


def savefig(filename, fig=None, scp=None, bbox_inches='tight', format=None, *args, **kwargs):
    assert format is None

    if filename[-4:].lower() not in ['.png', '.pdf', '.eps', '.svg'] and filename[-3:].lower() not in ['.ps']:
        extensions = ['.pdf', '.png', '.svg']
    else:
        extensions = ['']
    
    for ext in extensions:
        filename_with_ext = '{}{}'.format(filename, ext)
        if fig is None:
            import matplotlib.pylab as plt
            plt.savefig(filename_with_ext, bbox_inches=bbox_inches, format=None, *args, **kwargs)
        else:
            fig.savefig(filename_with_ext, bbox_inches=bbox_inches, format=None, *args, **kwargs)
        if scp is not None:
            scp_copy(source=filename_with_ext, destination=scp)


def inverse_transform(array_like, transform, offset):
    if transform is not None:
        return np.power(array_like, 1.0 / transform) - (0 if offset is None else offset)
    else:
        return array_like


def boxcox_transform(array_like):
    offset = ((-1 * np.min(array_like)) + 1) if np.min(array_like) <= 0 else None
    transformed_array, transform = boxcox(array_like + (0 if offset is None else offset))
    return transformed_array, transform, offset


def logit(scores):
    return np.exp(scores) / (1 + np.exp(scores))


def make_bold(text):
    return '$\mathcal{%s}$' % text


def make_italics(text):
    return '$\mathit{%s}$' % text


def make_bold_italics(text):
    return '$\mathbf{%s}$' % text


def sentence_case(s):
    return s[0].upper() + s[1:].lower()


def get_config(fn='config.yaml'):
    with open(fn, 'r') as f:
        config = yaml.load(f)
    return config


def simplify_dict(d, prefix=''):
    new_d = {}
    for k in d:
        if isinstance(d[k], dict):
            for kk in d[k]:
                new_d['{}{}_{}'.format(prefix, k, kk)] = d[k][kk]
        else:
            new_d[k] = d[k]
    return new_d


def is_unique(array_like):
    return np.unique(array_like).shape == np.array(array_like).flatten().shape


def no_nulls(array_like):
    return not np.any(pd.isnull(array_like))


def insert_str(source_str, insert_str, pos):
    return '{}{}{}'.format(source_str[:pos], insert_str, source_str[pos:])


def column_shape(a):
    if len(a.shape) == 1:
        a = a.reshape((a.shape[0], 1))
    return a


def is_iterable(x):
    return isinstance(x, collections.Iterable) and not isinstance(x, string_types)


def nans_like(arr):
    return np.ones_like(arr) * np.nan


def nans(shape):
    return np.ones(shape) * np.nan


def fill_in_axis(arr):
    return arr if len(arr.shape) > 1 else np.reshape(arr, (-1, 1))


def categorical_accuracy(true_y, pred_y):
    assert true_y.shape == pred_y.shape
    assert len(true_y.shape) == 2
    assert np.all(np.sum(true_y, axis=1) == 1)
    assert np.all(np.sum(true_y == 1, axis=1) == 1)
    return np.sum(np.argmax(true_y, axis=1) == np.argmax(pred_y, axis=1)) / float(len(true_y))


def mse(y_true, y_pred, axis=None):
    y_true = np.array(y_true)
    y_pred = np.array(y_pred)
    not_nan = ~pd.isnull(y_true)
    return np.mean(np.square(y_true[not_nan] - y_pred[not_nan]), axis=axis)


def cross_val_predict(estimator, X, y=None, groups=None, cv=None, n_jobs=1, verbose=0,
                      fit_params=None, pre_dispatch='2*n_jobs', method='predict', positive_class=1):
    pred = _cross_val_predict(estimator=estimator, X=X, y=y, groups=groups, cv=cv, n_jobs=n_jobs, verbose=verbose,
                              fit_params=fit_params, pre_dispatch=pre_dispatch, method=method)
    if method == 'predict_proba':
        pred = get_pred_proba(pred_proba=pred, classes=estimator.classes_, which_class=positive_class)
    return pred


def get_pred_proba(pred_proba, classes, which_class=1):
    return np.reshape(pred_proba[:, classes == which_class], (-1,))


def median_stratify(y, small=0, large=1):
    strat = np.ones_like(y) * small
    strat[y >= np.median(y)] = large
    return strat


def maybe_stratified_splits(X, y, n_splits, shuffle=True, random_state=0):
    n_smallest_class = np.min([np.sum(y == label) for label in np.unique(y)])
    kfold_method = StratifiedKFold if n_smallest_class >= n_splits else KFold
    kfold = kfold_method(n_splits=n_splits, shuffle=shuffle, random_state=random_state)
    splits = kfold.split(X=X, y=y)
    return list(splits)


def quant_norm_df(df, samples_in_rows=False):
    assert not samples_in_rows
    rank_mean = df.stack().groupby(df.rank(method='first').stack().astype(int)).mean()
    df_qn = df.rank(method='min').stack().astype(int).map(rank_mean).unstack()
    return df_qn


def anndata_to_df(adata, raw=False, index_name='cell', columns_name='gene',
                  assert_finite=True, assert_not_null=True, include_obs=True, include_var=True,
                  copy_X=True):
    df = adata.to_df() if not raw else pd.DataFrame(adata.raw.X, index=adata.obs_names, columns=adata.var_names)
    assert not assert_finite or np.all(np.isfinite(df.values))
    assert not assert_not_null or np.all(~pd.isnull(df.values))
    df.index.name = index_name
    df.columns.name = columns_name

    if include_obs:
        df.index = pd.MultiIndex.from_arrays(
            np.concatenate([df.index.values.reshape(1, -1), adata.obs.values.T]),
            names=[df.index.name] + adata.obs.columns.tolist())

    if include_var:
        df.columns = pd.MultiIndex.from_arrays(
            np.concatenate([df.columns.values.reshape(1, -1), adata.var.values.T]),
            names=[df.columns.name] + adata.var.columns.tolist())

    return df.copy() if copy_X else df


def tpm(counts_matrix, transcript_len=1, norm_factors=1, log=False, pseudocount=0.5, libsize_psuedocount=1, samples_as_rows=False, million=1e6):
    """
    Transcript per million normalization of gene expression
    If transcript_len=1 then you get counts per million
    Setting pseudocount=0.5, libsize_psuedocount=1 ensures results equivalent to LIMMA-voom
    """
    if samples_as_rows:
        counts_matrix = counts_matrix.T

    rpk = counts_matrix.astype(float) / transcript_len
    tpm = million * (rpk + pseudocount) / (rpk.sum(axis=0) * norm_factors + libsize_psuedocount)
    tpm = tpm.T if samples_as_rows else tpm

    return np.log2(tpm) if log else tpm


def filter_by_reads(counts_df, min_group_n, cpm_df=None, min_count=10, min_total_count=15, large_min_n=10,
                    large_min_n_proportion=0.7, for_large_min_n_use_only_proportion=False,
                    samples_as_rows=False, verbose=True):
    if samples_as_rows:
        counts_df = counts_df.T
        if cpm_df is not None:
            cpm_df = cpm_df.T

    if cpm_df is None:
        cpm_df = tpm(counts_df, transcript_len=1, norm_factors=1, log=False,
                     pseudocount=0.5, libsize_psuedocount=0, samples_as_rows=False)

    if min_group_n > large_min_n:
        if for_large_min_n_use_only_proportion:
            min_n = min_group_n * large_min_n_proportion
        else:
            min_n = large_min_n + (min_group_n - large_min_n) * large_min_n_proportion
    else:
        min_n = min_group_n

    _million, _tolerance = 1e6, 1e-14
    median_lib_size = counts_df.sum(axis=0).median()
    cpm_cutoff = min_count / median_lib_size * _million
    keep_cpm = (cpm_df >= cpm_cutoff).sum(axis=1) >= min_n - _tolerance
    keep_min_total_count = counts_df.sum(axis=1) >= min_total_count - _tolerance
    if verbose:
        print('for_large_min_n_use_only_proportion', for_large_min_n_use_only_proportion)
        print('min_n', min_n)
        print('median_lib_size', median_lib_size)
        print('cpm_cutoff', cpm_cutoff)
        print('remove based on cpm_cutoff', (~keep_cpm).sum())
        print('additionally remove based on keep_min_total_count', (~keep_min_total_count[keep_cpm]).sum())
    return keep_cpm & keep_min_total_count


def normalize_per_cell(counts_matrix, cells_as_rows=True, counts_per_cell_after=1e4):
    """
    Normalize total counts per cell.

    :param counts_matrix: "cells x genes" or "genes x cells" matrix, see cells_as_rows
    :param cells_as_rows: "True" if samples are as rows; "False" if samples are as columns
    :param counts_per_cell_after: after normalization, each cell will have a total count equal to this number
    :return: per-cell-normalized matrix
    """
    if cells_as_rows:
        counts_matrix = counts_matrix.T

    total_per_cell = counts_matrix.sum(axis=0).astype(float)
    normed_counts = counts_matrix / total_per_cell * counts_per_cell_after
    assert np.all(np.isnan(normed_counts.loc[:, total_per_cell == 0]))
    if np.any(total_per_cell == 0):
        print('WARNING: setting normalized counts of cells with no genes expressed to 0')
    normed_counts.loc[:, total_per_cell == 0] = 0

    return normed_counts.T if cells_as_rows else normed_counts


def collect_esat_counts(samples):
    """
    Collect gene expression (read counts, gene-level) output from ESAT into expression matrix for `samples`.

    samples : list of tuples (sample_name, sample_counts_fn)
    """

    expr_df = pd.DataFrame()
    for sample_name, sample_counts_fn in samples:
        with open_file(sample_counts_fn) as f:
            counts_df = pd.read_csv(f, sep='\t', usecols=['Symbol', 'Exp1'])
        counts_df = counts_df.rename(columns={'Symbol': 'gene_symbol', 'Exp1': sample_name}).set_index('gene_symbol')
        expr_df = pd.concat([expr_df, counts_df], axis=1, sort=True)

    return expr_df.sort_index(axis=1).fillna(0)


def get_palette(palette=None, name=None, n_colors=None, desat=None, start=None, rot=None, light=None, dark=None, reverse=False, as_cmap=False):
    if palette == 'color':
        if n_colors is None:
            n_colors = 10
        assert name in [None, 'deep', 'muted', 'bright', 'pastel', 'dark', 'colorblind']
        assert start is None and rot is None and light is None and dark is None and not reverse and not as_cmap
        SEABORN_PALETTE_LEN = 10
        rearrange_red = list(range(SEABORN_PALETTE_LEN))
        rearrange_red[1] = 3
        rearrange_red[3] = 1
        palette = sns.color_palette(palette=name, n_colors=SEABORN_PALETTE_LEN, desat=desat)
        pal_cycle = itertools.cycle(np.array(palette)[rearrange_red])
        palette = [next(pal_cycle) for _ in range(n_colors)]
    elif palette == 'sequential':
        assert n_colors is not None
        assert name is None and desat is None
        start = start if start is not None else 0
        rot = rot if rot is not None else 0.4
        light = light if light is not None else 0.85
        dark = dark if dark is not None else 0.15
        palette = sns.cubehelix_palette(n_colors=n_colors, start=start, rot=rot, light=light, dark=dark, reverse=reverse, as_cmap=as_cmap)
    elif palette == 'diverging':
        assert n_colors is not None
        assert start is None and rot is None and light is None and dark is None and not reverse
        assert name in [None, 'Spectral_r', 'Spectral', 'BrBG_r', 'BrBG', 'RdBu_r', 'RdBu']
        palette = sns.color_palette('Spectral_r' if name is None else name, n_colors=n_colors, desat=desat)
        if as_cmap:
            palette = ListedColormap(palette.as_hex())
    else:
        # e.g. palette='cubehelix'
        palette = sns.color_palette(palette=palette, n_colors=n_colors, desat=desat)
        if as_cmap:
            palette = ListedColormap(palette.as_hex())
    return palette


def clustermap_annot_colors(index, annotations, names=None, palettes=None, special_colors={}):
    assert names is None or len(annotations) == len(names)
    assert palettes is None or isinstance(palettes, string_types) or len(palettes) == len(annotations)
    assert np.all(len(index) == np.array([len(a) for a in annotations])), (len(index), np.array([len(a) for a in annotations]))

    if palettes is None:
        palettes = len(annotations) * ['color']
    elif isinstance(palettes, string_types):
        palettes = len(annotations) * [palettes]

    row_colors, color_codes = [], []
    for i, labels in enumerate(annotations):
        uniq_labels = sorted(set(labels).difference(special_colors))
        palette = get_palette(palette=palettes[i], n_colors=len(uniq_labels))
        color_codes.append(dict(zip(uniq_labels, palette)))
        for c in special_colors:
            if c in labels.astype(object):
                color_codes[-1][c] = special_colors[c]
        name = names[i] if names is not None else None
        # labels cannot be Series with an index, thus convert to array
        row_colors.append(pd.Series(np.array(labels), name=name, index=index).map(color_codes[-1]))
    return pd.concat(row_colors, axis=1), color_codes


def permute_one_axis_demo(a=None, idx=None):
    if a is None:
        a = np.random.randint(low=0, high=99, size=(100,7))

    if idx is None:
        idx = np.array([np.random.permutation(range(a.shape[1])) for i in range(a.shape[0])])

    assert a.shape == idx.shape
    add = np.arange(0, a.shape[0] * a.shape[1], a.shape[1])
    b = a.flatten()[((idx.T + add).T).flatten()].reshape(a.shape)

    assert np.all(a.sum(axis=1) == b.sum(axis=1))

    return a, idx, b


def make_dir(*dir_path):
    fig_dir = ''
    for d in dir_path:
        fig_dir = os.path.join(fig_dir, d)
        if not os.path.exists(fig_dir):
            try:
                os.mkdir(fig_dir)
            except FileExistsError:
                pass
    return fig_dir


def one_split(s, sep, last=False):
    s_split = s.split(sep)
    return s_split[0 if not last else -1], sep.join(s_split[1:] if not last else s_split[:-1])


def constant(a, axis=0):
    ref = a[0, :] if axis == 0 else a[:, 0][:, None]
    return np.all(a == ref, axis=axis)


def upper_triangle(df, k=1, other=np.nan):
    return df.where(np.triu(np.ones_like(df), k=k).astype(bool), other=other)


def lower_triangle(df, k=-1, other=np.nan):
    return df.where(np.tril(np.ones_like(df), k=k).astype(bool), other=other)


def non_diagonal(df):
    return df.mask(diagonal_mask(df))


def diagonal(df):
    return df.where(diagonal_mask(df))


def diagonal_mask(df):
    assert df.shape[0] == df.shape[1]
    return np.eye(df.shape[0], dtype=bool)


# this is equivalent to my VarianceSelector: SelectKBest(lambda _X, _y: _X.var(axis=0), n_components)
# class VarianceSelector(BaseEstimator, SelectorMixin):
#
#     def __init__(self, n_components):
#         self.n_components = n_components
#
#     def fit(self, X, y=None):
#         X = check_array(X)
#         self.variances_ = np.var(X, axis=0)
#         self.features_ = np.argsort(self.variances_)[::-1][:self.n_components]
#         return self
#
#     def _get_support_mask(self):
#         check_is_fitted(self, 'features_')
#         return np.in1d(range(len(self.variances_)), self.features_)


def argsort(seq):
    return sorted(range(len(seq)), key=seq.__getitem__)


def custom_argsort(seq, sort_dict):
    return sorted(range(len(seq)), key=lambda x: sort_dict[seq[x]])


def annotate_scatter(annotations, xy, ax):
    for ann, xy_pos in zip(annotations, xy):
        ax.annotate(ann, xy_pos)


def annotate_adata_scatter(annotation, adata, plot_type, ax):
    not_nan = ~pd.isnull(adata.obs[annotation])
    annotate_scatter(adata.obs[annotation].loc[not_nan], adata.obsm['X_{}'.format(plot_type)][not_nan, :2], ax)


def multirank_genes(pvals, LFCs, confs, pval_smaller_better=True, LFC_smaller_better=False, confs_smaller_better=False, group_ranks_func=np.max):
    assert np.min(pvals) >= 0 and np.max(pvals) <= 1
    m = np.array([np.array(pvals) if pval_smaller_better else -np.array(pvals),
                  np.array(LFCs) if LFC_smaller_better else -np.array(LFCs),
                  np.array(confs) if confs_smaller_better else -np.array(confs)])
    return multirank_matrix(m, group_ranks_func=group_ranks_func)


def multirank_matrix(m, smaller_better=True, group_ranks_func=np.max):
    if not smaller_better:
        m = - np.array(m)
    ranks = [rankdata(row, method='max') for row in m]
    return group_ranks_func(ranks, axis=0)


def plot_sequenced_samples(df, n_samples=30, ax=None, xlabel='Read count', ylabel='Density', title=None, xlim=None, ylim=None,
                           kde=True, hist=False, legend=False, samples_as_rows=False):
    if samples_as_rows:
        df = df.T
    for i in np.random.choice(list(range(df.shape[1])), size=n_samples, replace=False) if n_samples is not None else range(df.shape[1]):
        sample = df.iloc[:, i]
        ax = sns.distplot(sample, ax=ax, kde=kde, hist=hist, label=df.columns[i] if legend else None)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    if xlim is not None:
        ax.set_xlim(xlim)
    if ylim is not None:
        ax.set_ylim(ylim)
    if title is None:
        ax.set_title('{} samples{}{}'.format(n_samples, ', xlim: {}'.format(xlim) if xlim is not None else '', ', ylim: {}'.format(ylim) if ylim is not None else ''))
    else:
        ax.set_title(title)
    if legend:
        ax.legend(bbox_to_anchor=(1, 1))
    sns.despine()
    return ax


def feature_variance(X, y=None):
    return np.var(X, axis=0)


def _corrwith(X, y, method='pearson', mean=False):
    if y.ndim == 1:
        return pd.DataFrame(data=X).corrwith(pd.Series(y), method=method).values
    else:
        corrs = nans((X.shape[1], y.shape[1]))
        for i in range(y.shape[1]):
            corrs[:, i] = pd.DataFrame(data=X).corrwith(pd.Series(y[:, i]), method=method).values
        if mean:
            corrs = np.mean(corrs, axis=1)
        assert corrs.ndim == 1 and len(corrs) == X.shape[1]
        return corrs


def feature_correlation(X, y, method='pearson', absolute=True,
                        n_bootstrap=None, n_low_corr_variance=None, choice_size=0.5,
                        random_state=None):
    assert n_bootstrap is None or n_low_corr_variance is None
    if n_bootstrap:
        rhos = []
        rs = np.random.RandomState(random_state)
        for _ in range(n_bootstrap):
            idx = rs.choice(np.arange(len(X)), size=len(X), replace=True)
            rhos.append(_corrwith(X=X[idx], y=y[idx], method=method, mean=True))
        rho = np.mean(rhos, axis=0)
    elif n_low_corr_variance:
        rhos = []
        rs = np.random.RandomState(random_state)
        for _ in range(n_low_corr_variance):
            idx = rs.choice(np.arange(len(X)), size=int(np.ceil(len(X) * choice_size)), replace=False)
            rhos.append(_corrwith(X=X[idx], y=y[idx], method=method, mean=True))
        qrho = np.mean(rhos, axis=0) * (1 / np.var(rhos, axis=0))
    else:
        rho = _corrwith(X=X, y=y, method=method, mean=True)
    return np.abs(rho) if absolute else rho


def f_classif_with_low_variance(X, y, n_repeats=1000, fraction=0.5, q=90, p_thr=0.1, random_state=None):
    pvals = np.ones((n_repeats, X.shape[1])) * np.nan
    rs = np.random.RandomState(random_state)
    for i in range(n_repeats):
        positive_class = (y == 1)
        pos_idx = rs.choice(np.arange(len(X[positive_class])), size=int(len(X[positive_class]) * fraction), replace=False)
        neg_idx = rs.choice(np.arange(len(X[~positive_class])), size=int(len(X[~positive_class]) * fraction), replace=False)
        pvals[i] = f_classif(X=np.concatenate([X[positive_class][pos_idx], X[~positive_class][neg_idx]]),
                             y=np.concatenate([y[positive_class][pos_idx], y[~positive_class][neg_idx]]))[1]
    mask = (np.percentile(pvals, q, axis=0) < p_thr)
    pvals = -np.log10(np.median(pvals, axis=0))
    pvals[~mask] = 0
    return pvals


def stability_f1000_score(X, y, n_repeats=1000, fraction=0.5, top_n=100):
    scores = np.ones((n_repeats, X.shape[1])) * np.nan
    # rs = np.random.RandomState(random_state)
    for i in range(n_repeats):
        positive_class = (y == 1)
        pos_idx = np.random.choice(np.arange(len(X[positive_class])), size=int(len(X[positive_class]) * fraction), replace=False)
        neg_idx = np.random.choice(np.arange(len(X[~positive_class])), size=int(len(X[~positive_class]) * fraction), replace=False)
        score_func_ret = f_classif(X=np.concatenate([X[positive_class][pos_idx], X[~positive_class][neg_idx]]),
                                    y=np.concatenate([y[positive_class][pos_idx], y[~positive_class][neg_idx]]))
        scores[i] = score_func_ret[0] if isinstance(score_func_ret, (list, tuple)) else score_func_ret
    assert np.all(scores >= -1e-6)
    stability = (pd.DataFrame(data=scores).rank(axis=1, ascending=False) < top_n).sum(axis=0)
    return stability


def stability_f_score(X, y, n_repeats=100, fraction=0.5, top_n=100):
    scores = np.ones((n_repeats, X.shape[1])) * np.nan
    # rs = np.random.RandomState(random_state)
    for i in range(n_repeats):
        positive_class = (y == 1)
        pos_idx = np.random.choice(np.arange(len(X[positive_class])), size=int(len(X[positive_class]) * fraction), replace=False)
        neg_idx = np.random.choice(np.arange(len(X[~positive_class])), size=int(len(X[~positive_class]) * fraction), replace=False)
        score_func_ret = f_classif(X=np.concatenate([X[positive_class][pos_idx], X[~positive_class][neg_idx]]),
                                    y=np.concatenate([y[positive_class][pos_idx], y[~positive_class][neg_idx]]))
        scores[i] = score_func_ret[0] if isinstance(score_func_ret, (list, tuple)) else score_func_ret
    assert np.all(scores >= -1e-10)
    stability = (pd.DataFrame(data=scores).rank(axis=1, ascending=False) < top_n).sum(axis=0)
    return stability


def stability_selection(X, y, score_func, n_repeats=1000, fraction=0.5, top_n=100, random_state=None):
    scores = np.ones((n_repeats, X.shape[1])) * np.nan
    rs = np.random.RandomState(random_state)
    for i in range(n_repeats):
        positive_class = (y == 1)
        pos_idx = rs.choice(np.arange(len(X[positive_class])), size=int(len(X[positive_class]) * fraction), replace=False)
        neg_idx = rs.choice(np.arange(len(X[~positive_class])), size=int(len(X[~positive_class]) * fraction), replace=False)
        score_func_ret = score_func(X=np.concatenate([X[positive_class][pos_idx], X[~positive_class][neg_idx]]),
                                    y=np.concatenate([y[positive_class][pos_idx], y[~positive_class][neg_idx]]))
        scores[i] = score_func_ret[0] if isinstance(score_func_ret, (list, tuple)) else score_func_ret
    assert np.all(scores >= -1e-10)
    stability = (pd.DataFrame(data=scores).rank(axis=1, ascending=False) < top_n).sum(axis=0)
    return stability


def _get_de_peaks(model, k, seed, target):
    from misc import de_fn
    fn = de_fn(celltype='PBMC', model=model.format(target=target, seed=seed, k=k + 1))
    return pd.read_csv(fn, index_col=0)


def blind_DE_pval(contrast_suffix, fn=None, model=None, k=None, seed=None, target=None, peaks=None, region_filter=None):
    if fn:
        PEAK_ANNOT_ALL_FN = os.path.join('..', 'data', 'DE', 'peaks_filtered_PBMC.csv.gz')
        de_peaks_df = pd.read_csv(fn, index_col=0)
        if region_filter:
            peaks_df = pd.read_csv(PEAK_ANNOT_ALL_FN, index_col=0)
            assert de_peaks_df.index.equals(peaks_df.index)
            _peaks_mask = peaks_to_genes(peaks_df, **PEAKS_TO_GENES[region_filter])
            print(peaks_df.shape[0], _peaks_mask.sum())
            de_peaks_df.loc[~_peaks_mask, de_peaks_df.columns.str.contains('^p\.value\.')] = 1
        target = _TRIM_FN_PATTERN.search(fn)[1].replace('_FC1.2_responder', '')
    else:
        de_peaks_df = _get_de_peaks(model, k, seed, target).loc[peaks]
    return -np.log10(de_peaks_df.loc[:, 'p.value.{}{}'.format(target, contrast_suffix)].values)


def hashed_blind_DE_pval(X, contrast_suffix, target, region_filter=None):
    from misc import ML_RESULTS_ROOT
    HASH_DB = pd.read_csv(os.path.join(ML_RESULTS_ROOT, 'test_samples_nested_{}_{}.{}.hash_db.csv'.format(
        10, 10, target)), index_col=0)['blind_DE_fn']
    scores = blind_DE_pval(contrast_suffix=contrast_suffix, fn=HASH_DB[hash_sample_ids(X[:, 0])], region_filter=region_filter)
    scores = np.concatenate([[np.min(scores) - 1], scores], axis=0)
    assert np.argsort(scores)[0] == 0
    return scores


_TRIM_TARGET = 'thm.innate_nonspecific_24h_wo_LAC_IL10_IL1ra_V3'
_TRIM_FN_PATTERN = re.compile(r'de_results_p5_PBMC\.test_samples_nested_[0-9]+_[0-9]+\.([A-Za-z0-9\._]+)\.seed[0-9]+_.*\.csv(?:\.gz){0,1}$')


def hashed_blind_DE_N_V2_pval(X, y=None, target=_TRIM_TARGET):
    return hashed_blind_DE_pval(X=X, contrast_suffix='_FC1.2_N.V2', target=target)


def hashed_blind_DE_N_V3_pval(X, y=None, target=_TRIM_TARGET):
    return hashed_blind_DE_pval(X=X, contrast_suffix='_FC1.2_N.V3', target=target)


def hashed_blind_DE_R_V2_pval(X, y=None, target=_TRIM_TARGET):
    return hashed_blind_DE_pval(X=X, contrast_suffix='_FC1.2_R.V2', target=target)


def hashed_blind_DE_R_V3_pval(X, y=None, target=_TRIM_TARGET):
    return hashed_blind_DE_pval(X=X, contrast_suffix='_FC1.2_R.V3', target=target)


def hashed_blind_DE_N_pval(X, y=None, target=_TRIM_TARGET):
    return hashed_blind_DE_pval(X=X, contrast_suffix='_FC1.2_N', target=target)


def hashed_rand_DE_N_V2_pval(X, y=None, target=_TRIM_TARGET):
    return hashed_blind_DE_pval(X=X, contrast_suffix='_FC1.2_N.V2', target=target)


def hashed_rand_DE_N_V3_pval(X, y=None, target=_TRIM_TARGET):
    return hashed_blind_DE_pval(X=X, contrast_suffix='_FC1.2_N.V3', target=target)


def hashed_rand_DE_R_V2_pval(X, y=None, target=_TRIM_TARGET):
    return hashed_blind_DE_pval(X=X, contrast_suffix='_FC1.2_R.V2', target=target)


def hashed_rand_DE_R_V3_pval(X, y=None, target=_TRIM_TARGET):
    return hashed_blind_DE_pval(X=X, contrast_suffix='_FC1.2_R.V3', target=target)


def hashed_rand_DE_N_pval(X, y=None, target=_TRIM_TARGET):
    return hashed_blind_DE_pval(X=X, contrast_suffix='_FC1.2_N', target=target)


def hashed_shuffle_DE_N_V2_pval(X, y=None, target=_TRIM_TARGET):
    scores = hashed_blind_DE_pval(X=X, contrast_suffix='_FC1.2_N.V2', target=target)
    scores = np.random.RandomState(None).permutation(scores)
    scores[0] = np.min(scores) - 1
    return scores


def hashed_shuffle_DE_N_V3_pval(X, y=None, target=_TRIM_TARGET):
    scores = hashed_blind_DE_pval(X=X, contrast_suffix='_FC1.2_N.V3', target=target)
    scores = np.random.RandomState(None).permutation(scores)
    scores[0] = np.min(scores) - 1
    return scores


def hashed_shuffle_DE_R_V2_pval(X, y=None, target=_TRIM_TARGET):
    scores = hashed_blind_DE_pval(X=X, contrast_suffix='_FC1.2_R.V2', target=target)
    scores = np.random.RandomState(None).permutation(scores)
    scores[0] = np.min(scores) - 1
    return scores


def hashed_shuffle_DE_R_V3_pval(X, y=None, target=_TRIM_TARGET):
    scores = hashed_blind_DE_pval(X=X, contrast_suffix='_FC1.2_R.V3', target=target)
    scores = np.random.RandomState(None).permutation(scores)
    scores[0] = np.min(scores) - 1
    return scores


def hashed_shuffle_DE_N_pval(X, y=None, target=_TRIM_TARGET):
    scores = hashed_blind_DE_pval(X=X, contrast_suffix='_FC1.2_N', target=target)
    scores = np.random.RandomState(None).permutation(scores)
    scores[0] = np.min(scores) - 1
    return scores


def hashed_blind_DE_N_V2_pval_GENE_AND_DISTAL_10kb(X, y=None, target=_TRIM_TARGET):
    return hashed_blind_DE_pval(X=X, contrast_suffix='_FC1.2_N.V2', target=target, region_filter=GENE_AND_DISTAL_10kb)


def hashed_blind_DE_N_V3_pval_GENE_AND_DISTAL_10kb(X, y=None, target=_TRIM_TARGET):
    return hashed_blind_DE_pval(X=X, contrast_suffix='_FC1.2_N.V3', target=target, region_filter=GENE_AND_DISTAL_10kb)


def hashed_blind_DE_R_V2_pval_GENE_AND_DISTAL_10kb(X, y=None, target=_TRIM_TARGET):
    return hashed_blind_DE_pval(X=X, contrast_suffix='_FC1.2_R.V2', target=target, region_filter=GENE_AND_DISTAL_10kb)


def hashed_blind_DE_R_V3_pval_GENE_AND_DISTAL_10kb(X, y=None, target=_TRIM_TARGET):
    return hashed_blind_DE_pval(X=X, contrast_suffix='_FC1.2_R.V3', target=target, region_filter=GENE_AND_DISTAL_10kb)


def hashed_blind_DE_N_pval_GENE_AND_DISTAL_10kb(X, y=None, target=_TRIM_TARGET):
    return hashed_blind_DE_pval(X=X, contrast_suffix='_FC1.2_N', target=target, region_filter=GENE_AND_DISTAL_10kb)


def random_features(X, y=None):
    scores = np.random.rand(X.shape[1] - 1)
    scores = np.concatenate([[np.min(scores) - 1], scores], axis=0)
    assert np.argsort(scores)[0] == 0
    return scores


PICKABLE_SCORES = {
    'mutual_info_classif': mutual_info_classif,
    'mutual_info_regression': mutual_info_regression,
    'f_classif': f_classif,
    'f_regression': f_regression,
    'variance': feature_variance,
    'f_classif_stability': stability_f_score,
    'f1000_classif_stability': stability_f1000_score,

    'random': random_features,

    'blind_DE_N_V2_pval': hashed_blind_DE_N_V2_pval,
    'blind_DE_N_V3_pval': hashed_blind_DE_N_V3_pval,
    'blind_DE_R_V2_pval': hashed_blind_DE_R_V2_pval,
    'blind_DE_R_V3_pval': hashed_blind_DE_R_V3_pval,
    'blind_DE_N_pval':  hashed_blind_DE_N_pval,
    
    'rand_DE_N_V2_pval_1': hashed_rand_DE_N_V2_pval,
    'rand_DE_N_V3_pval_1': hashed_rand_DE_N_V3_pval,
    'rand_DE_R_V2_pval_1': hashed_rand_DE_R_V2_pval,
    'rand_DE_R_V3_pval_1': hashed_rand_DE_R_V3_pval,
    'rand_DE_N_pval_1': hashed_rand_DE_N_pval,
    
    'rand_DE_N_V2_pval_2': hashed_rand_DE_N_V2_pval,
    'rand_DE_N_V3_pval_2': hashed_rand_DE_N_V3_pval,
    'rand_DE_R_V2_pval_2': hashed_rand_DE_R_V2_pval,
    'rand_DE_R_V3_pval_2': hashed_rand_DE_R_V3_pval,
    'rand_DE_N_pval_2': hashed_rand_DE_N_pval,
    
    'rand_DE_N_V2_pval_3': hashed_rand_DE_N_V2_pval,
    'rand_DE_N_V3_pval_3': hashed_rand_DE_N_V3_pval,
    'rand_DE_R_V2_pval_3': hashed_rand_DE_R_V2_pval,
    'rand_DE_R_V3_pval_3': hashed_rand_DE_R_V3_pval,
    'rand_DE_N_pval_3': hashed_rand_DE_N_pval,
    
    'rand_DE_N_V2_pval_4': hashed_rand_DE_N_V2_pval,
    'rand_DE_N_V3_pval_4': hashed_rand_DE_N_V3_pval,
    'rand_DE_R_V2_pval_4': hashed_rand_DE_R_V2_pval,
    'rand_DE_R_V3_pval_4': hashed_rand_DE_R_V3_pval,
    'rand_DE_N_pval_4': hashed_rand_DE_N_pval,
    
    'rand_DE_N_V2_pval_5': hashed_rand_DE_N_V2_pval,
    'rand_DE_N_V3_pval_5': hashed_rand_DE_N_V3_pval,
    'rand_DE_R_V2_pval_5': hashed_rand_DE_R_V2_pval,
    'rand_DE_R_V3_pval_5': hashed_rand_DE_R_V3_pval,
    'rand_DE_N_pval_5': hashed_rand_DE_N_pval,
    
    'rand_DE_N_V2_pval_6': hashed_rand_DE_N_V2_pval,
    'rand_DE_N_V3_pval_6': hashed_rand_DE_N_V3_pval,
    'rand_DE_R_V2_pval_6': hashed_rand_DE_R_V2_pval,
    'rand_DE_R_V3_pval_6': hashed_rand_DE_R_V3_pval,
    'rand_DE_N_pval_6': hashed_rand_DE_N_pval,
    
    'rand_DE_N_V2_pval_7': hashed_rand_DE_N_V2_pval,
    'rand_DE_N_V3_pval_7': hashed_rand_DE_N_V3_pval,
    'rand_DE_R_V2_pval_7': hashed_rand_DE_R_V2_pval,
    'rand_DE_R_V3_pval_7': hashed_rand_DE_R_V3_pval,
    'rand_DE_N_pval_7': hashed_rand_DE_N_pval,
    
    'rand_DE_N_V2_pval_8': hashed_rand_DE_N_V2_pval,
    'rand_DE_N_V3_pval_8': hashed_rand_DE_N_V3_pval,
    'rand_DE_R_V2_pval_8': hashed_rand_DE_R_V2_pval,
    'rand_DE_R_V3_pval_8': hashed_rand_DE_R_V3_pval,
    'rand_DE_N_pval_8': hashed_rand_DE_N_pval,
    
    'rand_DE_N_V2_pval_9': hashed_rand_DE_N_V2_pval,
    'rand_DE_N_V3_pval_9': hashed_rand_DE_N_V3_pval,
    'rand_DE_R_V2_pval_9': hashed_rand_DE_R_V2_pval,
    'rand_DE_R_V3_pval_9': hashed_rand_DE_R_V3_pval,
    'rand_DE_N_pval_9': hashed_rand_DE_N_pval,
    
    'rand_DE_N_V2_pval_10': hashed_rand_DE_N_V2_pval,
    'rand_DE_N_V3_pval_10': hashed_rand_DE_N_V3_pval,
    'rand_DE_R_V2_pval_10': hashed_rand_DE_R_V2_pval,
    'rand_DE_R_V3_pval_10': hashed_rand_DE_R_V3_pval,
    'rand_DE_N_pval_10': hashed_rand_DE_N_pval,
    
    'shuffle_rand_DE_N_V2_pval_1': hashed_shuffle_DE_N_V2_pval,
    'shuffle_rand_DE_N_V3_pval_1': hashed_shuffle_DE_N_V3_pval,
    'shuffle_rand_DE_R_V2_pval_1': hashed_shuffle_DE_R_V2_pval,
    'shuffle_rand_DE_R_V3_pval_1': hashed_shuffle_DE_R_V3_pval,
    'shuffle_rand_DE_N_pval_1': hashed_shuffle_DE_N_pval,
    
    'shuffle_rand_DE_N_V2_pval_2': hashed_shuffle_DE_N_V2_pval,
    'shuffle_rand_DE_N_V3_pval_2': hashed_shuffle_DE_N_V3_pval,
    'shuffle_rand_DE_R_V2_pval_2': hashed_shuffle_DE_R_V2_pval,
    'shuffle_rand_DE_R_V3_pval_2': hashed_shuffle_DE_R_V3_pval,
    'shuffle_rand_DE_N_pval_2': hashed_shuffle_DE_N_pval,
    
    'shuffle_rand_DE_N_V2_pval_3': hashed_shuffle_DE_N_V2_pval,
    'shuffle_rand_DE_N_V3_pval_3': hashed_shuffle_DE_N_V3_pval,
    'shuffle_rand_DE_R_V2_pval_3': hashed_shuffle_DE_R_V2_pval,
    'shuffle_rand_DE_R_V3_pval_3': hashed_shuffle_DE_R_V3_pval,
    'shuffle_rand_DE_N_pval_3': hashed_shuffle_DE_N_pval,
    
    'shuffle_rand_DE_N_V2_pval_4': hashed_shuffle_DE_N_V2_pval,
    'shuffle_rand_DE_N_V3_pval_4': hashed_shuffle_DE_N_V3_pval,
    'shuffle_rand_DE_R_V2_pval_4': hashed_shuffle_DE_R_V2_pval,
    'shuffle_rand_DE_R_V3_pval_4': hashed_shuffle_DE_R_V3_pval,
    'shuffle_rand_DE_N_pval_4': hashed_shuffle_DE_N_pval,
    
    'shuffle_rand_DE_N_V2_pval_5': hashed_shuffle_DE_N_V2_pval,
    'shuffle_rand_DE_N_V3_pval_5': hashed_shuffle_DE_N_V3_pval,
    'shuffle_rand_DE_R_V2_pval_5': hashed_shuffle_DE_R_V2_pval,
    'shuffle_rand_DE_R_V3_pval_5': hashed_shuffle_DE_R_V3_pval,
    'shuffle_rand_DE_N_pval_5': hashed_shuffle_DE_N_pval,
    
    'shuffle_rand_DE_N_V2_pval_6': hashed_shuffle_DE_N_V2_pval,
    'shuffle_rand_DE_N_V3_pval_6': hashed_shuffle_DE_N_V3_pval,
    'shuffle_rand_DE_R_V2_pval_6': hashed_shuffle_DE_R_V2_pval,
    'shuffle_rand_DE_R_V3_pval_6': hashed_shuffle_DE_R_V3_pval,
    'shuffle_rand_DE_N_pval_6': hashed_shuffle_DE_N_pval,
    
    'shuffle_rand_DE_N_V2_pval_7': hashed_shuffle_DE_N_V2_pval,
    'shuffle_rand_DE_N_V3_pval_7': hashed_shuffle_DE_N_V3_pval,
    'shuffle_rand_DE_R_V2_pval_7': hashed_shuffle_DE_R_V2_pval,
    'shuffle_rand_DE_R_V3_pval_7': hashed_shuffle_DE_R_V3_pval,
    'shuffle_rand_DE_N_pval_7': hashed_shuffle_DE_N_pval,
    
    'shuffle_rand_DE_N_V2_pval_8': hashed_shuffle_DE_N_V2_pval,
    'shuffle_rand_DE_N_V3_pval_8': hashed_shuffle_DE_N_V3_pval,
    'shuffle_rand_DE_R_V2_pval_8': hashed_shuffle_DE_R_V2_pval,
    'shuffle_rand_DE_R_V3_pval_8': hashed_shuffle_DE_R_V3_pval,
    'shuffle_rand_DE_N_pval_8': hashed_shuffle_DE_N_pval,
    
    'shuffle_rand_DE_N_V2_pval_9': hashed_shuffle_DE_N_V2_pval,
    'shuffle_rand_DE_N_V3_pval_9': hashed_shuffle_DE_N_V3_pval,
    'shuffle_rand_DE_R_V2_pval_9': hashed_shuffle_DE_R_V2_pval,
    'shuffle_rand_DE_R_V3_pval_9': hashed_shuffle_DE_R_V3_pval,
    'shuffle_rand_DE_N_pval_9': hashed_shuffle_DE_N_pval,
    
    'shuffle_rand_DE_N_V2_pval_10': hashed_shuffle_DE_N_V2_pval,
    'shuffle_rand_DE_N_V3_pval_10': hashed_shuffle_DE_N_V3_pval,
    'shuffle_rand_DE_R_V2_pval_10': hashed_shuffle_DE_R_V2_pval,
    'shuffle_rand_DE_R_V3_pval_10': hashed_shuffle_DE_R_V3_pval,
    'shuffle_rand_DE_N_pval_10': hashed_shuffle_DE_N_pval,
    
    'blind_DE_N_V2_pval_GENE_AND_DISTAL_10kb': hashed_blind_DE_N_V2_pval_GENE_AND_DISTAL_10kb,
    'blind_DE_N_V3_pval_GENE_AND_DISTAL_10kb': hashed_blind_DE_N_V3_pval_GENE_AND_DISTAL_10kb,
    'blind_DE_R_V2_pval_GENE_AND_DISTAL_10kb': hashed_blind_DE_R_V2_pval_GENE_AND_DISTAL_10kb,
    'blind_DE_R_V3_pval_GENE_AND_DISTAL_10kb': hashed_blind_DE_R_V3_pval_GENE_AND_DISTAL_10kb,
    'blind_DE_N_pval_GENE_AND_DISTAL_10kb':  hashed_blind_DE_N_pval_GENE_AND_DISTAL_10kb
}


FEATURE_SELECTION_SCORES = {
    'mutual_info_classif': mutual_info_classif,
    'mutual_info_regression': mutual_info_regression,
    'f_classif': f_classif,
    'f_regression': f_regression,
    'variance': feature_variance,
    'spearman_rho': lambda X, y: feature_correlation(
        X, y, method='spearman', absolute=True, n_bootstrap=None),
    'spearman_rho_with_bootstrap': lambda X, y, seed: feature_correlation(
        X, y, method='spearman', absolute=True, n_bootstrap=10, random_state=seed),
    'f_classif_with_low_variance': lambda X, y, p_thr, seed: f_classif_with_low_variance(
        X, y, n_repeats=1000, fraction=0.5, q=90, p_thr=p_thr, random_state=seed),
    'f_classif_stability': lambda X, y, top_n, seed: stability_selection(
        X, y, score_func=f_classif, n_repeats=100, fraction=0.5, top_n=top_n, random_state=seed),
    # 'f_classif_100_stability': lambda X, y, top_n, seed: stability_selection(
    #     X, y, score_func=f_classif, n_repeats=100, fraction=0.5, top_n=top_n, random_state=seed),
    'f_regression_stability': lambda X, y, top_n, seed: stability_selection(
        X, y, score_func=f_regression, n_repeats=100, fraction=0.5, top_n=top_n, random_state=seed),
    'blind_DE_N_V2_pval': lambda k, seed, target, peaks: blind_DE_pval(contrast_suffix='_N.V2', model='test_samples.{target}.seed{seed}_k{k}', k=k, seed=seed, target=target, peaks=peaks),
    'blind_DE_N_V3_pval': lambda k, seed, target, peaks: blind_DE_pval(contrast_suffix='_N.V3', model='test_samples.{target}.seed{seed}_k{k}', k=k, seed=seed, target=target, peaks=peaks),
    'blind_DE_R_V2_pval': lambda k, seed, target, peaks: blind_DE_pval(contrast_suffix='_R.V2', model='test_samples.{target}.seed{seed}_k{k}', k=k, seed=seed, target=target, peaks=peaks),
    'blind_DE_R_V3_pval': lambda k, seed, target, peaks: blind_DE_pval(contrast_suffix='_R.V3', model='test_samples.{target}.seed{seed}_k{k}', k=k, seed=seed, target=target, peaks=peaks),
    'blind_DE_N_pval': lambda k, seed, target, peaks: blind_DE_pval(contrast_suffix='_N', model='test_samples.{target}.seed{seed}_k{k}', k=k, seed=seed, target=target, peaks=peaks),
}


def binarize_labels(f, binary_mapper=None):
    """ Decorator useful for binarizing labels (y variable).
    It works well with fit and score methods of estimators
    as well as with __call__ method of scorers.
    """
    @wraps(f)
    def wrapper(estimator, X, y, *args, **kwargs):
        assert isinstance(estimator, BaseEstimator)
        if y is not None:
            assert X.shape[0] == y.shape[0]
            y = binary_mapper(y) if binary_mapper is not None else estimator.binary_mapper(y)
        return f(estimator, X, y, *args, **kwargs)

    return wrapper


class BinarizedLogisticRegression(LogisticRegression):

    def __init__(self, binary_mapper,
                 penalty='l2', dual=False, tol=1e-4, C=1.0,
                 fit_intercept=True, intercept_scaling=1, class_weight=None,
                 random_state=None, solver='warn', max_iter=100,
                 multi_class='warn', verbose=0, warm_start=False, n_jobs=None,
                 l1_ratio=None):
        super(BinarizedLogisticRegression, self).__init__(
            penalty=penalty, dual=dual, tol=tol, C=C, fit_intercept=fit_intercept, intercept_scaling=intercept_scaling,
            class_weight=class_weight, random_state=random_state, solver=solver,
            max_iter=max_iter, multi_class=multi_class, verbose=verbose, warm_start=warm_start, n_jobs=n_jobs,
            l1_ratio=l1_ratio)
        self.binary_mapper = binary_mapper

    @binarize_labels
    def fit(self, X, y, *args, **kwargs):
        super(BinarizedLogisticRegression, self).fit(X, y, *args, **kwargs)

    @binarize_labels
    def score(self, X, y, *args, **kwargs):
        super(BinarizedLogisticRegression, self).score(X, y, *args, **kwargs)


class BinarizedRandomForestClassifier(RandomForestClassifier):

    def __init__(self, binary_mapper,
                 n_estimators='warn',
                 criterion="gini",
                 max_depth=None,
                 min_samples_split=2,
                 min_samples_leaf=1,
                 min_weight_fraction_leaf=0.,
                 max_features="auto",
                 max_leaf_nodes=None,
                 min_impurity_decrease=0.,
                 min_impurity_split=None,
                 bootstrap=True,
                 oob_score=False,
                 n_jobs=None,
                 random_state=None,
                 verbose=0,
                 warm_start=False,
                 class_weight=None
                 ):
        super(BinarizedRandomForestClassifier, self).__init__(
            n_estimators=n_estimators,
            criterion=criterion,
            max_depth=max_depth,
            min_samples_split=min_samples_split,
            min_samples_leaf=min_samples_leaf,
            min_weight_fraction_leaf=min_weight_fraction_leaf,
            max_features=max_features,
            max_leaf_nodes=max_leaf_nodes,
            min_impurity_decrease=min_impurity_decrease,
            min_impurity_split=min_impurity_split,
            bootstrap=bootstrap,
            oob_score=oob_score,
            n_jobs=n_jobs,
            random_state=random_state,
            verbose=verbose,
            warm_start=warm_start,
            class_weight=class_weight)
        self.binary_mapper = binary_mapper

    @binarize_labels
    def fit(self, X, y, *args, **kwargs):
        super(BinarizedRandomForestClassifier, self).fit(X, y, *args, **kwargs)

    @binarize_labels
    def score(self, X, y, *args, **kwargs):
        super(BinarizedRandomForestClassifier, self).score(X, y, *args, **kwargs)


def _read_list_dict_or_none(obj, index, key, default=None):
    if obj is None:
        return None
    elif isinstance(obj, dict):
        assert key is not None
        return obj.get(key, default)
    else:
        return obj[index]


def dim_reduced_scatterplots(df, n_dims=2, labels=None, palettes=None, legends='brief',
                             orders=None, hue_orders='sort', num_legend_fmt='auto',
                             xlabel=None, ylabel=None, title=None, figsize=(6, 3), xlim=None, ylim=None,
                             wspace=0.8, hspace=0.5, sharex=False, sharey=False, show_xticks=True, show_yticks=True,
                             xticks=None, yticks=None, axes_style='white',
                             legend_transformer=None, legend_markerscale=None, legend_ncol=1, legend_fontsize=None,
                             random_state=0, **scatter_kwargs):
    """
    :param df: typically a pca_df, must have a MultiIndex with names being the labels to plot
    :param n_dims: how many top dimensions (PCs) to plot
    :param labels: subset of the names in df's MultiIndex to plot
    :param palettes: string or list or dict of palettes for the labels
    :param legends: string or list or dict of legends for the labels
    :param orders: string or list or dict of orders for the labels
    :param hue_orders: string or list or dict of hue_orders for the labels
    :param num_legend_fmt: a way to correct number format in legend (use if 'auto' does not work well)
    :param xlabel: xlabel
    :param ylabel: ylabel
    :param title: title
    :param figsize: figsize
    :param xlim: xlim
    :param ylim: ylim
    :param wspace: the amount of width reserved for space between subplots, expressed as a fraction of the average axis width
    :param hspace: the amount of height reserved for space between subplots, expressed as a fraction of the average axis height
    :param sharex: bool or {'none', 'all', 'row', 'col'}, controls sharing of properties among x axes
    :param sharey: bool or {'none', 'all', 'row', 'col'}, controls sharing of properties among y axes
    :param axes_style: seaborn style: 'darkgrid', 'whitegrid', 'dark', 'white' or 'ticks'
    :param random_state: random_state
    :param scatter_kwargs: anything matplotlib.pyplot.scatter accepts
    :return: matplotlib.figure.Figure and array of matplotlib.axes.Axes objects tuple
    """

    assert labels is not None or (palettes is None and orders is None)
    if labels is not None and isinstance(labels, str):
        labels = [labels]

    if palettes is not None and isinstance(palettes, str):
        palettes = [palettes] * (1 if labels is None else len(labels))
    assert palettes is None or len(palettes) == len(labels) or isinstance(palettes, dict)

    if orders is not None and isinstance(orders, str):
        orders = [orders] * (1 if labels is None else len(labels))
    assert orders is None or len(orders) == len(labels) or isinstance(orders, dict)

    if hue_orders is None or isinstance(hue_orders, str):
        hue_orders = [hue_orders] * (1 if labels is None else len(labels))
    assert (labels is None and len(hue_orders) == 1) or (
                labels is not None and (len(hue_orders) == len(labels) or isinstance(hue_orders, dict)))

    if legends is None or isinstance(legends, str):
        legends = [legends] * (1 if labels is None else len(labels))
    assert (labels is None and len(legends) == 1) or (
                labels is not None and (len(legends) == len(labels) or isinstance(legends, dict)))

    nrows = n_dims - 1
    ncols = 1 if labels is None else len(labels)

    with sns.axes_style(axes_style):
        fig, axs = plt.subplots(nrows=nrows, ncols=ncols, squeeze=False, sharex=sharex, sharey=sharey,
                                figsize=(figsize[0] * ncols, figsize[1] * nrows))
        plt.subplots_adjust(wspace=wspace, hspace=hspace)
        for row in range(axs.shape[0]):
            for col in range(axs.shape[1]):
                ax = axs[row, col]
                _label = labels[col] if labels is not None else None
                _palette = _read_list_dict_or_none(palettes, col, _label)
                _legend = _read_list_dict_or_none(legends, col, _label, default='brief')
                _order = _read_list_dict_or_none(orders, col, _label)

                if _order is not None:
                    if _order == 'random':
                        np.random.seed(random_state)
                        df = df.iloc[np.random.permutation(df.shape[0])]
                    else:
                        assert _order in ['ascending', 'descending']
                        assert _label is not None
                        df = df.sort_index(level=_label, ascending=(_order == 'ascending'))

                _hue = df.index.get_level_values(_label).values if _label is not None else None
                _hue_order = _read_list_dict_or_none(hue_orders, col, _label,
                                                     default='sort') if _label is not None else None

                sns.scatterplot(ax=ax,
                                x=df.iloc[:, row].values,
                                y=df.iloc[:, row + 1].values,
                                hue=_hue,
                                hue_order=np.unique(_hue) if _hue_order == 'sort' else _hue_order,
                                palette=_palette,
                                legend=_legend,
                                **scatter_kwargs)
                sns.despine()

                if _label is not None and _legend is not None:
                    _handles, _labels = ax.get_legend_handles_labels()
                    if legend_transformer:
                        print('Legend before transforming:', _labels)
                        _labels = [legend_transformer(l) for l in _labels]
                    legend = ax.legend(_handles, _labels, bbox_to_anchor=(1, 1), title=_label,
                                       markerscale=legend_markerscale, ncol=legend_ncol, fontsize=legend_fontsize)
                    if num_legend_fmt != 'auto' and _hue.dtype == float:
                        for t in legend.texts:
                            t.set_text(('{:' + num_legend_fmt + '}').format(float(t.get_text())))
                if xlim is not None:
                    ax.set_xlim(xlim)
                if ylim is not None:
                    ax.set_ylim(ylim)
                ax.set_xlabel('{} {}'.format(xlabel if xlabel is not None else 'Dimension', row + 1))
                ax.set_ylabel('{} {}'.format(ylabel if ylabel is not None else 'Dimension', row + 2))
                if not show_xticks:
                    ax.set_xticklabels([])
                if not show_yticks:
                    ax.set_yticklabels([])
                if xticks is not None:
                    ax.set_xticks(xticks)
                if yticks is not None:
                    ax.set_yticks(yticks)
                ax.set_title(_label)
        if title is not None:
            fig.suptitle(title)

    return fig, axs


def explained_var_plot(pca, expl_var_thr):
    n_components_above_thr = np.sum(pca.explained_variance_ratio_ > expl_var_thr) if expl_var_thr <= 1 else expl_var_thr
    ax = sns.pointplot(x=list(range(1, n_components_above_thr + 1)),
                       y=pca.explained_variance_ratio_[:n_components_above_thr])
    ax.set_xlabel('Principal component ID')
    ax.set_ylabel('Explained variance ratio')
    ax.set_title('Principal components with explained\nvariance ratio > {}'.format(expl_var_thr))
    sns.despine()
    return ax


def force_df_to_numeric(df, include_binary=True, verbose=True):
    assert isinstance(df, pd.DataFrame)
    df = df.copy()
    for col in df:
        try:
            df[col] = pd.to_numeric(df[col])
        except ValueError as e:
            if verbose:
                print('{} is not numeric: {}'.format(col, e))
            remove = not include_binary
            if include_binary:
                df[col] = df[col].astype('category')
                if len(df[col].cat.categories) == 2:
                    try:
                        cat_map = {}
                        for cat in df[col].cat.categories:
                            cat_map[cat] = df[col].cat.codes[df[col] == cat][0]
                    except Exception as e:
                        cat_map = str(e)
                    assert df[col].isnull().equals(df[col].cat.codes == -1)
                    if verbose:
                        print('Converting {} to categorical: {}'.format(col, cat_map))
                    df[col] = df[col].cat.codes
                    df[col] = df[col].where(df[col] != -1).astype(float)
                else:
                    if verbose:
                        print('{} has more than two categories: {}{}'.format(
                        col, ', '.join(df[col].cat.categories[:4]), ', ...' if len(df[col].cat.categories) > 4 else ''))
                    remove = True

            if remove:
                if verbose:
                    print('Removing {}'.format(col))
                df = df.drop(col, axis=1)
    return df


def keep_categorical(df, verbose=True):
    assert isinstance(df, pd.DataFrame)
    for col in df:
        if not hasattr(df[col].dtype, 'cat'):
            try:
                pd.to_numeric(df[col])
            except ValueError as e:
                pass
            else:
                # Remove those that can be converted to numeric
                if verbose:
                    print('Removing {}'.format(col))
                df = df.drop(col, axis=1)
    return df.astype('category') if df.shape[1] != 0 else df


def test_association(df, categories_df, axis=0, test='nonparametric', return_test_stat=False, verbose=True):
    assert isinstance(df, pd.DataFrame) and isinstance(categories_df, pd.DataFrame)
    assert not isinstance(test, str) or test in ['nonparametric', 'parametric']
    if axis == 1:
        df, categories_df = df.T, categories_df.T
    assert df.index.equals(categories_df.index)

    skipped = []
    assoc_df = pd.DataFrame(index=categories_df.columns, columns=df.columns)
    for col in df.columns:
        values = df[col].values
        for category_name in categories_df.columns:
            labels = categories_df[category_name].values
            not_nan = ~pd.isnull(values) & ~pd.isnull(labels)
            categories = np.unique(labels[not_nan])
            if len(categories) > 1:
                if isinstance(test, str):
                    if len(categories) == 2:
                        _test_method = ttest_ind if test == 'parametric' else mannwhitneyu
                    else:
                        _test_method = f_oneway if test == 'parametric' else kruskal
                else:
                    _test_method = test
                try:
                    test_stat, pval = _test_method(*[values[not_nan][labels[not_nan] == c] for c in categories])
                except DeprecationWarning as w:
                    print(w)
                    pass
                except ValueError:
                    result = np.nan
                else:
                    result = pval if not return_test_stat else test_stat
            else:
                if category_name not in skipped:
                    if verbose:
                        print('Skipping {} because it has too few categories: {}'.format(category_name, categories))
                    skipped.append(category_name)
                result = np.nan
            assoc_df.loc[category_name, col] = result

    return assoc_df.astype(float)


def test_contingency(df1, df2, axis=0, test=chi2_contingency, return_test_stat=False, verbose=True):
    assert isinstance(df1, pd.DataFrame) and isinstance(df2, pd.DataFrame)
    if axis == 1:
        df1, df2 = df1.T, df2.T
    assert df1.index.equals(df2.index)

    skipped = []
    assoc_df = pd.DataFrame(index=df2.columns, columns=df1.columns)
    for col1 in df1.columns:
        values1 = df1[col1].astype('category')
        for col2 in df2.columns:
            values2 = df2[col2].astype('category')
            if len(values1.cat.categories) > 1 and len(values2.cat.categories) > 1:
                result = test(pd.crosstab(values1.reset_index(drop=True), values2.reset_index(drop=True)))
                result = result[1] if not return_test_stat else result[0]
            else:
                for col, n_cat in [(col1, len(values1.cat.categories)), (col2, len(values2.cat.categories))]:
                    if n_cat <= 1 and col not in skipped:
                        if verbose:
                            print('Skipping {} because it has too few categories'.format(col))
                        skipped.append(col)
                result = np.nan
            assoc_df.loc[col2, col1] = result

    return assoc_df.astype(float)


def correlate_dataframes(df1, df2, axis=0):
    assert isinstance(df1, pd.DataFrame) and isinstance(df2, pd.DataFrame)
    if axis == 1:
        df1, df2 = df1.T, df2.T
    assert df1.index.equals(df2.index)

    corr_df = pd.DataFrame(columns=df1.columns)
    for col in df2.columns:
        corr_df = corr_df.append(df1.corrwith(df2[col], axis=0).rename(col))
    return corr_df


def dim_reduced_assoc_heatmap(df=None, n_dims=None, labels=None, assoc_df=None,
                              data_type='numeric', annot_type='numeric', binary_as_numeric=True, assoc_test='nonparametric',
                              assoc_test_stat=False,
                              annot_as_rows=True, abs_corr=False, verbose=True, **heatmap_kws):
    """
    :param df: typically a pca_df, must have a MultiIndex with names being the labels to plot
    :param n_dims: how many top dimensions (PCs) to plot
    :param labels: subset of the names in df's MultiIndex to plot
    :param assoc_df: if df is None, use precalculated association scores
    :param annot_type: depending on the selected labels 'numeric' or 'category'
    :param binary_as_numeric: include binary categorical labels as numeric
    :param assoc_test: 'parametric' or 'nonparametric' or method(*arrays)
    :param assoc_test_stat: if annot_type is 'category', plot test statistics instead of p-values
    :param annot_as_rows: annotations displayed as rows, data dimensions (PCs) as columns
    :param abs_corr: if annot_type is 'numeric', plot absolute correlation coefficients
    :param heatmap_kws: anything seaborn.heatmap can accept
    :return: matplotlib.axes.Axes object
    """

    assert 'data' not in heatmap_kws
    assert data_type in ['numeric', 'category']
    assert annot_type in ['numeric', 'category']
    if assoc_df is None:
        annot_df = df.index.to_frame()
        if n_dims is not None:
            df = df.iloc[:, :n_dims]
        if labels is not None:
            annot_df = annot_df.loc[:, labels]

        if data_type == 'numeric' and annot_type == 'numeric':
            assoc_df = correlate_dataframes(df, force_df_to_numeric(annot_df, include_binary=binary_as_numeric, verbose=verbose))
        elif data_type != annot_type:
            assoc_df = test_association(df, keep_categorical(annot_df, verbose=verbose), test=assoc_test, return_test_stat=assoc_test_stat, verbose=verbose)
        else:
            assert data_type == 'category' and annot_type == 'category'
            assoc_df = test_contingency(keep_categorical(df, verbose=verbose), keep_categorical(annot_df, verbose=verbose), return_test_stat=assoc_test_stat, verbose=verbose)

        if not annot_as_rows:
            assoc_df = assoc_df.T
    else:
        if n_dims is not None:
            assoc_df = assoc_df.iloc[:, :n_dims] if annot_as_rows else assoc_df.iloc[:n_dims]
        if labels is not None:
            assoc_df = assoc_df.loc[labels] if annot_as_rows else assoc_df.loc[:, labels]

    if abs_corr:
        assoc_df = assoc_df.abs()

    if 'ax' not in heatmap_kws:
        heatmap_kws['ax'] = plt.subplots(nrows=1, ncols=1,
                                         figsize=(0.9 * len(assoc_df.columns), 0.9 * len(assoc_df.index)))[1]
    if 'cmap' not in heatmap_kws and annot_type == 'category' and not assoc_test_stat:
        heatmap_kws['cmap'] = 'rocket_r'

    ax = sns.heatmap(assoc_df, **heatmap_kws)
    ax.set_ylim(len(assoc_df), 0)
    ax.set_yticklabels(ax.get_yticklabels(), rotation=0)
    if data_type == 'numeric' and annot_type == 'numeric':
        ax.set_title('Pearson correlation')
    else:
        ax.set_title('{}{}'.format('P-value' if not assoc_test_stat else 'Test statistic',
                                   ' ({})'.format(assoc_test) if data_type != annot_type else ''))

    return ax


def pca_analysis(categorical_labels=None, continuous_labels=None, df=None, pca_df=None, annot_df=None, pca=None, scale=True,
                 force_title=None, title_prefix=None, fig_dir='', fig_prefix='', fig_ext='svg', save_fig=True, dpi=300, random_state=0,
                 expl_var_thr=0.01,
                 plot_variance_explained=True, plot_scatter=True, verbose=True, plot_assoc_heatmap=True,
                 try_categorical_with_correlation=True, show=False,
                 scatter_kwargs={}, categorical_heatmap_kwargs={}, continuous_heatmap_kwargs={}):
    assert categorical_labels or continuous_labels
    assert np.sum([df is not None, pca_df is not None]) == 1
    if categorical_labels is None:
        categorical_labels = {}
    if continuous_labels is None:
        continuous_labels = {}
    if isinstance(categorical_labels, str):
        categorical_labels = {categorical_labels: categorical_labels}
    elif not isinstance(categorical_labels, dict):
        categorical_labels = {None: categorical_labels}
    if isinstance(continuous_labels, str):
        continuous_labels = {continuous_labels: continuous_labels}
    elif not isinstance(continuous_labels, dict):
        continuous_labels = {None: continuous_labels}

    # performa PCA
    if pca_df is None:
        pca = PCA(random_state=random_state)
        pca_df = pd.DataFrame(data=pca.fit_transform(StandardScaler().fit_transform(df) if scale else df),
                              index=df.index)
        pca_df.columns = ['PC{}'.format(i + 1) for i in range(pca_df.shape[1])]
        categorical_pca_df = None
        assert pca_df.index.equals(df.index)
    else:
        pca_df = pca_df.copy()
        categorical_pca_df = keep_categorical(pca_df, verbose=False)

    # if annot_df is not None, create a MultiIndex from it
    if annot_df is not None:
        assert pca_df.index.equals(annot_df.index)
        pca_df.index = pd.MultiIndex.from_frame(annot_df)

    expl_var_ax, scatter_axs, heatmap_axs = None, None, None

    # variance explained
    if plot_variance_explained and pca is not None:
        expl_var_ax = explained_var_plot(pca, expl_var_thr=expl_var_thr)
        for _thr in [0.1, 0.01, 0.005, 0.001]:
            if _thr > expl_var_thr:
                expl_var_ax.axhline(_thr, ls='--', c='gray')
        if title_prefix is not None:
            expl_var_ax.set_title('{}: {}'.format(title_prefix, expl_var_ax.get_title()))
        if show or save_fig:
            if save_fig:
                savefig(os.path.join(fig_dir, '{}explained_var.{}'.format(fig_prefix, fig_ext)), dpi=dpi)
            if show:
                plt.show()
            plt.close()
    else:
        expl_var_ax = None

    # scatter plots
    if plot_scatter:
        scatter_axs = []
        for labels in [categorical_labels, continuous_labels]:
            for label_type in labels.keys():
                not_nan = np.ones((pca_df.shape[0],), dtype=bool)
                for c in labels[label_type]:
                    not_nan &= ~pca_df.index.get_level_values(c).isnull()
                if verbose:
                    print('Missing {} data:'.format(label_type), (~not_nan).sum())

                if 'random_state' not in scatter_kwargs:
                    scatter_kwargs['random_state'] = random_state
                fig, axs = dim_reduced_scatterplots(
                    pca_df[not_nan], labels=labels[label_type], **scatter_kwargs)
                if force_title:
                    for ax in axs.flatten():
                        ax.set_title(force_title)
                elif title_prefix is not None:
                    for ax in axs.flatten():
                        ax.set_title('{}: {}'.format(title_prefix, ax.get_title()))
                scatter_axs.append(axs)
                if show or save_fig:
                    if save_fig:
                        savefig(os.path.join(fig_dir, '{}scatter{}.{}'.format(fig_prefix, '.{}'.format(label_type) if label_type is not None else '', fig_ext)), dpi=dpi)
                    if show:
                        plt.show()
                    plt.close()

    # association heatmaps
    if plot_assoc_heatmap:
        _cat_columns, _cat_pca_tuple = (categorical_pca_df.columns, [('category', categorical_pca_df)]) \
            if categorical_pca_df is not None and categorical_pca_df.shape[1] != 0 else ([], [])
        for data_type, _pca_df in [('numeric', pca_df.drop(_cat_columns, axis=1))] + _cat_pca_tuple:
            categorical_heatmap_kwargs['annot_type'] = 'category'
            continuous_heatmap_kwargs['annot_type'] = 'numeric'
            heatmap_axs = []
            for labels, heatmap_kwargs in [
                (categorical_labels, categorical_heatmap_kwargs),
                # if pca_df is numeric, treat categorical labels with correlation as well
                (categorical_labels if try_categorical_with_correlation and data_type == 'numeric' else {}, continuous_heatmap_kwargs),
                (continuous_labels, continuous_heatmap_kwargs)
            ]:
                for label_type in labels.keys():
                    ax = dim_reduced_assoc_heatmap(
                        _pca_df, labels=labels[label_type], data_type=data_type, verbose=verbose, **heatmap_kwargs)
                    if title_prefix is not None:
                        ax.set_title('{}: {}'.format(title_prefix, ax.get_title()))
                    heatmap_axs.append(ax)
                    if save_fig or show:
                        if save_fig:
                            savefig(os.path.join(fig_dir, '{}heatmap.{}{}{}.{}'.format(
                                fig_prefix, heatmap_kwargs['annot_type'], '_vs_category' if data_type == 'category' else '',
                                '.{}'.format(label_type) if label_type is not None else '', fig_ext)), dpi=dpi)
                        if show:
                            plt.show()
                        plt.close()

    return pca_df, expl_var_ax, scatter_axs, heatmap_axs


def clustermap(data,
               row_anns=None, row_anns_palettes=None, row_anns_legend=None,
               col_anns=None, col_anns_palettes=None, col_anns_legend=None,
               pivot_kws=None, method='average', metric='euclidean',
               z_score=None, standard_scale=None, figsize=None, cbar_kws=None,
               row_cluster=True, col_cluster=True, row_linkage=None, col_linkage=None,
               row_colors=None, col_colors=None, mask=None,
               vmin=None, vmax=None, cmap=None, center=None, robust=False,
               xticklabels=False, yticklabels=False,
               xticklabels_kws=None, yticklabels_kws=None,
               default_palettes=None,
               show_row_dendrogram=True, show_col_dendrogram=True,
               row_legend_bbox=None, col_legend_bbox=None,
               row_legend_ncol=1, col_legend_ncol=1,
               row_anns_center=0, col_anns_center=0,
               hide_ticks=True, xlabel=None, ylabel=None,
               title=None, figure_fn=None, raster_heatmap=None, raster_row_dendrogram=None, raster_col_dendrogram=None,
               raster_row_colors=None, raster_col_colors=None, raster_cbar=None, **kwargs):

    assert z_score is None or z_score in [0, 1]

    for colors, anns, palettes, legend in [(row_colors, row_anns, row_anns_palettes, row_anns_legend),
                                           (col_colors, col_anns, col_anns_palettes, col_anns_legend)]:
        assert not (colors is not None and anns is not None)
        assert not (anns is None and (palettes is not None or legend is not None))
        assert not (anns is not None and palettes is not None and len(anns) != len(palettes))
        assert not (anns is not None and legend is not None and not set(legend).issubset(anns))

    if isinstance(data, pd.DataFrame):
        df = data
    elif isinstance(data, Iterable):
        # A list/tuple of dataframes
        assert len(data) in [1, 2, 3]
        df = data[0]
        if len(data) > 1:
            obs_df = data[1]
        if len(data) > 2:
            var_df = data[2]
    else:
        # AnnData or similar
        df = data.to_df()
        obs_df = data.obs
        var_df = data.var
    assert np.all(np.isfinite(df.values)) and np.all(~pd.isnull(df.values))

    keep_rows = ~constant(df.values, axis=1) if df.shape[1] > 1 else np.ones((df.shape[0],), dtype=bool)
    keep_cols = ~constant(df.loc[keep_rows].values, axis=0)

    if np.any(~keep_rows) or np.any(~keep_cols):
        print(
            'Removing {} rows and {} columns because they are constant'.format(np.sum(~keep_rows), np.sum(~keep_cols)))
        df = df.loc[keep_rows, keep_cols]
        if obs_df is not None:
            obs_df = obs_df.loc[obs_df.index[keep_rows]]
        if var_df is not None:
            var_df = var_df.loc[var_df.index[keep_cols]]

    if row_anns is not None:
        if row_anns_palettes is None and default_palettes is not None:
            row_anns_palettes = [default_palettes.get(a) for a in row_anns]
        row_colors, row_cc = clustermap_annot_colors(index=df.index, annotations=[obs_df[a] for a in row_anns],
                                                     ann_names=row_anns, palettes=row_anns_palettes,
                                                     center=row_anns_center)
    if col_anns is not None:
        if col_anns_palettes is None and default_palettes is not None:
            col_anns_palettes = [default_palettes.get(a) for a in col_anns]
        col_colors, col_cc = clustermap_annot_colors(index=df.columns, annotations=[var_df[a] for a in col_anns],
                                                     ann_names=col_anns, palettes=col_anns_palettes,
                                                     center=col_anns_center)

    cg = sns.clustermap(df,
                        pivot_kws=pivot_kws, method=method, metric=metric,
                        z_score=z_score, standard_scale=standard_scale, figsize=figsize, cbar_kws=cbar_kws,
                        row_cluster=row_cluster, col_cluster=col_cluster,
                        row_linkage=row_linkage, col_linkage=col_linkage,
                        row_colors=row_colors, col_colors=col_colors, mask=mask,
                        vmin=vmin, vmax=vmax, cmap=cmap, center=center, robust=robust,
                        xticklabels=xticklabels, yticklabels=yticklabels,
                        rasterized=False, **kwargs)
    cg.ax_heatmap.set_ylim(len(df), 0)

    if ylabel is not None:
        cg.ax_heatmap.set_ylabel(ylabel.format(n=df.shape[0]) if '{n}' in ylabel else ylabel)
    if xlabel is not None:
        cg.ax_heatmap.set_xlabel(xlabel.format(n=df.shape[1]) if '{n}' in xlabel else xlabel)
    if raster_heatmap:
        for heatmap_child in cg.ax_heatmap.get_children():
            if isinstance(heatmap_child, mpl.collections.QuadMesh):
                heatmap_child.set_rasterized(True)
    cg.cax.set_rasterized(raster_cbar)
    cg.ax_row_dendrogram.set_rasterized(raster_row_dendrogram)
    cg.ax_col_dendrogram.set_rasterized(raster_col_dendrogram)
    if cg.ax_row_colors is not None:
        cg.ax_row_colors.set_rasterized(raster_row_colors)
        # cg.ax_row_colors.set_ylim(len(row_anns), 0)
    if cg.ax_col_colors is not None:
        cg.ax_col_colors.set_rasterized(raster_col_colors)
        # cg.ax_col_colors.set_ylim(len(col_anns), 0)

    if title is not None:
        cg.fig.suptitle(title)
    cg.ax_row_dendrogram.set_visible(show_row_dendrogram)
    cg.ax_col_dendrogram.set_visible(show_col_dendrogram)
    if hide_ticks:
        cg.ax_heatmap.tick_params(axis='both', which='both', length=0)
    if xticklabels and xticklabels_kws is not None:
        cg.ax_heatmap.set_xticklabels(cg.ax_heatmap.get_xticklabels(), **xticklabels_kws)
    if yticklabels and yticklabels_kws is not None:
        cg.ax_heatmap.set_yticklabels(cg.ax_heatmap.get_yticklabels(), **yticklabels_kws)
    # cg.cax.set_position([-0.2, 0.15, 0.3, 0.03])  # [x, y, width, height]

    if row_anns is not None:
        cg.ax_row_colors.set_xticklabels(cg.ax_row_colors.get_xticklabels(), fontsize=max(2, 11 - len(row_anns)))
        cg.ax_row_colors.tick_params(axis='both', which='both', length=0)

    if col_anns is not None:
        cg.ax_col_colors.set_yticklabels(cg.ax_col_colors.get_yticklabels(), fontsize=max(2, 11 - len(col_anns)))
        cg.ax_col_colors.tick_params(axis='both', which='both', length=0)

    merged_legend_titles = not ((row_anns_legend is not None and len(row_anns_legend) > 1) or \
                           (col_anns_legend is not None and len(col_anns_legend) > 1))

    if row_anns_legend is not None:
        row_anns = np.array(row_anns)
        row_anns_palettes = np.array(row_anns_palettes)
        row_legend_palettes = [row_anns_palettes[row_anns == a][0] for a in row_anns_legend]
        legend_for_clustermap_annot_colors(cg.ax_row_dendrogram, ann_names=row_anns_legend, color_codes=row_cc,
                                           palettes=row_legend_palettes, bbox_to_anchor=row_legend_bbox,
                                           ncol=row_legend_ncol, merged_legend_titles=merged_legend_titles,
                                           center=row_anns_center)

    if col_anns_legend is not None:
        col_anns = np.array(col_anns)
        col_anns_palettes = np.array(col_anns_palettes)
        col_legend_palettes = [col_anns_palettes[col_anns == a][0] for a in col_anns_legend]
        legend_for_clustermap_annot_colors(cg.ax_col_dendrogram, ann_names=col_anns_legend, color_codes=col_cc,
                                           palettes=col_legend_palettes, bbox_to_anchor=col_legend_bbox,
                                           ncol=col_legend_ncol, merged_legend_titles=merged_legend_titles,
                                           center=col_anns_center)

    if figure_fn is not None:
        savefig(figure_fn)

    return cg


def clustermap_annot_colors(index, annotations, ann_names, palettes=None, center=0, special_colors={}):
    '''
    :param index: pd.Dataframe.index for which the colors to be plotted
    :param annotations: 2-D array-like values for different annotations
    :param ann_names: 1-D array-like names of the annotations
    :param palettes: None, str, or list of str; 'S' for the default sequential, 'D' for the default 'diverging',
    and 'Q' for the default qualitative palette. Or a specific palette name accepted by sns.color_palette (note
    that in this case, data will be treated in a categorical/qualitative fashion)
    :param special_colors: dictionary to replace e.g. np.nan with default color
    :return: pd.Dataframe with colors for sns.clustermap, and dict of dicts with color codes for every annotation
    '''

    S = ['S', 'sequential']
    D = ['D', 'diverging']

    assert len(annotations) == len(ann_names)
    assert palettes is None or isinstance(palettes, string_types) or len(palettes) == len(annotations)
    assert np.all(len(index) == np.array([len(a) for a in annotations])), (len(index), np.array([len(a) for a in annotations]))

    if palettes is None:
        palettes = len(annotations) * [None]
    elif isinstance(palettes, string_types):
        palettes = len(annotations) * [palettes]

    row_colors, color_codes = [], {}
    for i, (ann_name, values) in enumerate(zip(ann_names, annotations)):
        values = np.array(values)
        if pd.isnull(values).any():
            print('Warning: {} contains NaNs'.format(ann_name))

        if palettes[i] not in ['S', 'sequential', 'D', 'diverging']:
            uniq_values = np.unique(values)
            _palette = sns.color_palette(
                palette=palettes[i] if is_iterable(palettes[i]) else PALETTE_CODES.get(palettes[i], palettes[i]),
                n_colors=len(uniq_values))
            color_codes[ann_name] = dict(zip(uniq_values, _palette))

        else:
            _palette_name = PALETTE_CODES[palettes[i]]
            _cmap = plt.get_cmap(_palette_name) if _palette_name != 'cubehelix' else sns.cubehelix_palette(as_cmap=True)
            if isinstance(palettes[i], str) and palettes[i] in S:
                _vmin = values.min()
                _vmax = values.max()
            elif isinstance(palettes[i], str) and palettes[i] in D:
                vrange = max(values.max() - center, center - values.min())
                _vmin = center - vrange
                _vmax = center + vrange
            else:
                raise ValueError
            _norm = mpl.colors.Normalize(vmin=_vmin, vmax=_vmax)
            _colors = _cmap(_norm(values))
            color_codes[ann_name] = dict(zip(values, _colors))

        _replace_special_colors(color_codes[ann_name], special_colors)
        # values cannot be Series with an index, must be an array!
        row_colors.append(pd.Series(values, name=ann_name, index=index).map(color_codes[ann_name]))

    return pd.concat(row_colors, axis=1), color_codes


def _replace_special_colors(color_dict, special_colors):
    for value in special_colors:
        if value in color_dict:
            color_dict[value] = special_colors[value]


def legend_for_clustermap_annot_colors(dendrogram_ax, ann_names, color_codes, palettes=None, center=0,
                                       loc='best', bbox_to_anchor=None, ncol=1, adjust_col_length=True,
                                       merged_legend_titles=True, **kwargs):
    assert palettes is None or len(palettes) == len(ann_names)

    S = ['S', 'sequential']
    D = ['D', 'diverging']
    S_cls = [0, 25, 50, 75, 100]
    D_cls = [0, 50, 100]

    col_length = None
    if adjust_col_length and len(ann_names) > 1 and ncol > 1:
        col_length = np.max(
            [len(color_codes[ann_name]) if palettes is None or palettes[i] not in S + D else (
                len(S_cls) if isinstance(palettes[i], str) and palettes[i] in S else 2 * len(D_cls) - 1
            ) for i, ann_name in enumerate(ann_names)]
        )

    for i, ann_name in enumerate(ann_names):
        values = np.sort(list(color_codes[ann_name].keys()))
        if palettes is not None:
            if isinstance(palettes[i], str) and palettes[i] in S:
                values = np.percentile(values, S_cls, interpolation='nearest')
            elif isinstance(palettes[i], str) and palettes[i] in D:
                assert np.sum(values < center) > 1 and np.sum(values >= center) > 1
                values = np.concatenate([np.percentile(values[values < center], D_cls, interpolation='nearest'),
                                         np.percentile(values[values >= center], D_cls, interpolation='nearest')])
                assert sorted(values) == values.tolist() and len(values) == 6
                neg_is_closer_to_zero = abs(values[2]) < abs(values[3])
                values = values[[True, True, neg_is_closer_to_zero, not neg_is_closer_to_zero, True, True]]
        # if there are several annotations, insert empty lines between them (and possibly make titles)
        if len(ann_names) > 1 or not merged_legend_titles:
            dendrogram_ax.bar(0, 0, linewidth=0, color='white', label=ann_name if not merged_legend_titles else ' ')
        for value in values:
            dendrogram_ax.bar(0, 0, linewidth=0, color=color_codes[ann_name][value], label=value)
        if col_length is not None:
            for _ in range(col_length - len(values)):
                dendrogram_ax.bar(0, 0, linewidth=0, color='white', label=' ')
    dendrogram_ax.legend(loc=loc, bbox_to_anchor=bbox_to_anchor, ncol=ncol,
                         title=', '.join(ann_names) if merged_legend_titles else None,
                         **kwargs)
    # else:
    #     values = list(color_codes[ann_names[0]].keys())
    #     sns.scatterplot(x=pd.Series(values, index=values),
    #                     y=pd.Series(values, index=values),
    #                     hue=pd.Series(values, index=values),
    #                     palette=color_codes[ann_names[0]],
    #                     ax=dendrogram_ax)
    #     dendrogram_ax.legend(loc='best', bbox_to_anchor=[0, 1], frameon=True, title=ann_names[0])


def binarize(df, negative=(0.0, 0.25), positive=(0.75, 1.0), return_q=False):
    assert isinstance(df, pd.Series) or isinstance(df, pd.DataFrame)
    quantiles = df.quantile(q=sorted(set([negative[0], negative[1], positive[0], positive[1]])),
                               **(dict(axis=0) if isinstance(df, pd.DataFrame) else dict()))
    binary_df = pd.DataFrame(index=df.index, columns=df.columns, dtype=float) if isinstance(df, pd.DataFrame) else pd.Series(index=df.index, dtype=float)
    binary_df = binary_df.mask((df >= quantiles.loc[negative[0]]) & (df <= quantiles.loc[negative[1]]), other=0)
    binary_df = binary_df.mask((df >= quantiles.loc[positive[0]]) & (df <= quantiles.loc[positive[1]]), other=1)
    return binary_df if not return_q else (binary_df, quantiles)


def get_binarize_thresholds(df, negative=(0.0, 0.25), positive=(0.75, 1.0)):
    assert isinstance(df, pd.Series) or isinstance(df, pd.DataFrame)
    quantiles = df.quantile(q=sorted(set([negative[0], negative[1], positive[0], positive[1]])),
                               **(dict(axis=0) if isinstance(df, pd.DataFrame) else dict()))
    return quantiles


def binarize_with_thresholds(df, quantiles, negative=(0.0, 0.25), positive=(0.75, 1.0)):
    binary_df = pd.DataFrame(index=df.index, columns=df.columns, dtype=float) if isinstance(df, pd.DataFrame) else pd.Series(index=df.index, dtype=float)
    binary_df = binary_df.mask((df >= quantiles.loc[negative[0]]) & (df <= quantiles.loc[negative[1]]), other=0)
    binary_df = binary_df.mask((df >= quantiles.loc[positive[0]]) & (df <= quantiles.loc[positive[1]]), other=1)
    return binary_df


def mask_with_thresholds(df, quantiles, negative=(0.0, 0.25), positive=(0.75, 1.0)):
    masked_df = df.copy().astype(float)
    masked_df = masked_df.where(((df >= quantiles.loc[negative[0]]) & (df <= quantiles.loc[negative[1]])) | ((df >= quantiles.loc[positive[0]]) & (df <= quantiles.loc[positive[1]])))
    return masked_df


def _binarize_deprecated(df, negative=(0.0, 0.25), positive=(0.75, 1.0), return_q=False):
    assert isinstance(df, pd.Series) or isinstance(df, pd.DataFrame)
    quantiles = pd.DataFrame()
    for qs, interpolation in [([negative[0], positive[0]], 'lower'), ([negative[1], positive[1]], 'higher')]:
        quantiles = pd.concat([quantiles,
                               df.quantile(q=qs, interpolation=interpolation,
                                           **(dict(axis=0) if isinstance(df, pd.DataFrame) else dict()))], axis=0)
    quantiles = quantiles.loc[[negative[0], negative[1], positive[0], positive[1]]]
    if isinstance(df, pd.Series):
        assert quantiles.shape[1] == 1
        quantiles = quantiles[0].rename(df.name)
    binary_df = pd.DataFrame(index=df.index, columns=df.columns, dtype=float) if isinstance(df, pd.DataFrame) else pd.Series(index=df.index, dtype=float)
    binary_df = binary_df.mask((df >= quantiles.loc[negative[0]]) & (df <= quantiles.loc[negative[1]]), other=0)
    binary_df = binary_df.mask((df >= quantiles.loc[positive[0]]) & (df <= quantiles.loc[positive[1]]), other=1)
    return binary_df if not return_q else (binary_df, quantiles)


def tukeys_outliers(df, k=1.5):
    q = df.quantile([0.25, 0.75])
    lower = q.loc[0.25] - (k * (q.loc[0.75] - q.loc[0.25]))
    upper = q.loc[0.75] + (k * (q.loc[0.75] - q.loc[0.25]))
    return lower, upper


def _rename_param_grid_params(param_grid, prefix):
    if isinstance(param_grid, dict):
        return {'{}{}'.format(ESTIMATOR_PREFIX, p): param_grid[p] for p in param_grid}
    elif is_iterable(param_grid):
        return [{'{}{}'.format(ESTIMATOR_PREFIX, p): d[p] for p in d} for d in param_grid]


class SelectRandomKBest(SelectKBest):

    def __init__(self, score_func, k, variances, means, variances_rtol=0.1, means_rtol=0.1, remove_features=None, random_state=None):
        super().__init__(score_func=score_func, k=k)
        self.means = means
        self.variances = variances
        self.means_rtol = means_rtol
        self.variances_rtol = variances_rtol
        self.remove_features = remove_features
        self.random_state = random_state

    def _get_support_mask(self):
        check_is_fitted(self)

        if self.k == 'all':
            return np.ones(self.scores_.shape, dtype=bool)
        elif self.k == 0:
            return np.zeros(self.scores_.shape, dtype=bool)
        else:
            assert not np.isnan(self.scores_).any()

            already_selected = list(self.remove_features) if self.remove_features else []
            shuffled_scores = np.array(self.scores_)
            for idx in np.argsort(self.scores_, kind="mergesort")[-self.k:]:
                choose_from, i = [], 0
                while (len(choose_from) < 100) and (i < 10):
                    i += 1
                    variance_mask = (self.variances < self.variances[idx] * (1 + self.variances_rtol * i)) & (self.variances >= self.variances[idx] * (1 - self.variances_rtol * i))
                    mean_mask = (self.means < self.means[idx] * (1 + self.means_rtol * i)) & (self.means >= self.means[idx] * (1 - self.means_rtol * i))
                    available_mask = ~np.isin(np.arange(len(self.scores_)), already_selected)
                    choose_from = np.arange(len(self.scores_))[variance_mask & mean_mask & available_mask]
                if len(choose_from) < 100:
                    choose_from = np.arange(len(self.scores_))[available_mask]
                new_idx = np.random.RandomState(self.random_state).choice(choose_from, size=1)[0]
                score = self.scores_[idx]
                shuffled_scores[idx] = np.min(self.scores_)
                shuffled_scores[new_idx] = score
                already_selected.append(new_idx)
            self.scores_ = shuffled_scores
            sys.stdout.flush()
            mask = np.zeros(self.scores_.shape, dtype=bool)
            mask[np.argsort(self.scores_, kind="mergesort")[-self.k:]] = 1
            return mask


def grid_cv(estimator, X, y, groups, cv, param_grid, scorer, scale=True, feature_selection_score=None,
            n_features=None, random_features=False, n_pcs=None, pca_first=False, merged_folds_scoring=False,
            grid_train_scores=False, refit=True, n_jobs=None, verbose=1):

    if isinstance(estimator, RidgeClassifierCV):
        raise ValueError
        assert not grid_train_scores
        assert n_features is None
        assert n_pcs is None
        assert groups is None

        estimator.fit(X, y)

        best_score, best_alpha = None, None
        for i, alpha in enumerate(estimator.alphas):
            alpha_score = roc_auc_score(y, estimator.cv_values_[:, 0, i])
            if (best_score is None) or (alpha_score > best_score):
                best_score, best_alpha = alpha_score, alpha
        _all_alphas = list(estimator.alphas)
        estimator.set_params(alphas=[best_alpha])
        estimator.fit(X, y)
        estimator.set_params(alphas=_all_alphas)

        grid_search = GridSearchCV(estimator=estimator, param_grid=param_grid, scoring=scorer, refit=True, cv=cv,
                                   return_train_score=grid_train_scores, n_jobs=n_jobs)
        grid_search.best_estimator_ = estimator
        grid_search.best_score_ = best_score
        assert best_alpha == estimator.alpha_
        grid_search.best_params_ = {'alpha': best_alpha}

    elif isinstance(estimator, BaggingClassifier):
        raise ValueError
        print('This will be a BaggingClassifier')
        assert not grid_train_scores
        assert n_features is None
        assert n_pcs is None
        assert groups is None

        assert len(param_grid) == 1 and len(param_grid['n_estimators']) == 1, param_grid
        estimator.base_estimator.set_params(scoring=scorer)
        estimator.set_params(n_estimators=param_grid['n_estimators'][0], n_jobs=n_jobs)
        print(estimator)
        estimator.fit(X, y)

        grid_search = GridSearchCV(estimator=estimator, param_grid=param_grid, scoring=scorer, refit=True, cv=cv,
                                   return_train_score=grid_train_scores, n_jobs=n_jobs)
        grid_search.best_estimator_ = estimator
        grid_search.best_score_ = 0
        grid_search.best_params_ = {'n_estimators': estimator.n_estimators}

    else:
        steps, memory = [], False
        if pca_first:
            if scale:
                steps.append((SCALE_STEP, StandardScaler()))
            if n_pcs is not None:
                print('n_pcs', n_pcs if n_pcs > 0 else None)
                steps.append((PCA_STEP, PCA(n_components=n_pcs if n_pcs > 0 else None)))
                memory = True
        if n_features is not None:
            if random_features:
                # this is only true if applied to blind DE scores
                # because the first feature is in fact index in the hash database
                # I am making sure this is indeed the case
                fake_feature_idx = 0
                assert np.argsort(feature_selection_score(X, y))[0] == fake_feature_idx
                k_best_selector = SelectRandomKBest(score_func=feature_selection_score, k=n_features,
                                                    variances=X.var(axis=0), means=X.mean(axis=0),
                                                    variances_rtol=0.1, means_rtol=0.1,
                                                    remove_features=[fake_feature_idx])
            else:
                k_best_selector = SelectKBest(score_func=feature_selection_score, k=n_features)
            steps.append((FS_STEP, k_best_selector))
            memory = True
        if not pca_first:
            if scale:
                steps.append((SCALE_STEP, StandardScaler()))
            if n_pcs is not None:
                print('n_pcs', n_pcs if n_pcs > 0 else None)
                steps.append((PCA_STEP, PCA(n_components=n_pcs if n_pcs > 0 else None)))
                memory = True
        steps.append((EST_STEP, estimator))
        pipeline = Pipeline(steps, memory=joblib.Memory('./cache', verbose=verbose > 2) if memory else None,
                            verbose=verbose > 1)

        # merged_folds_scoring or LOO
        if merged_folds_scoring or len(cv) == len(X):
            assert is_iterable(cv)
            assert not grid_train_scores
            print('Grid search with merged_folds_scoring')
            best_score, best_params = None, None
            for params in ParameterGrid(_rename_param_grid_params(param_grid, ESTIMATOR_PREFIX)):
                pipeline.set_params(**params)
                y_pred, y_test = [], []
                for train, test in cv:
                    pipeline.fit(X[train], y[train])
                    if hasattr(pipeline, 'decision_function'):
                        y_pred.extend(pipeline.decision_function(X[test]))
                    elif hasattr(pipeline, 'predict_proba'):
                        y_pred.extend(pipeline.predict_proba(X[test])[:, pipeline.classes_ == 1].flatten())
                    else:
                        y_pred.extend(pipeline.predict(X[test]))
                    y_test.extend(y[test])
                assert len(y_pred) == len(y) and len(y_test) == len(y), (len(y_pred), len(y), len(y_test))
                score = get_scorer_func(scorer)(y_test, y_pred)
                if (best_score is None) or (score > best_score):
                    best_score, best_params = score, params

            # refit
            pipeline.set_params(**best_params)
            pipeline.fit(X, y)

            grid_search = GridSearchCV(estimator=pipeline,
                                       param_grid=_rename_param_grid_params(param_grid, ESTIMATOR_PREFIX),
                                       scoring=scorer, refit=True, cv=cv,
                                       return_train_score=grid_train_scores, n_jobs=n_jobs)
            grid_search.best_estimator_ = pipeline
            grid_search.best_score_ = best_score
            grid_search.best_params_ = best_params

        # classic CV (mean of folds' scores)
        else:
            print('Classic grid search')
            grid_search = GridSearchCV(estimator=pipeline,
                                       param_grid=_rename_param_grid_params(param_grid, ESTIMATOR_PREFIX),
                                       scoring=scorer, refit=refit, cv=cv,
                                       return_train_score=grid_train_scores, n_jobs=n_jobs)
            grid_search.fit(X, y, groups=groups)

    return grid_search


def make_integer_samples_ids(samples):
    return np.asarray([int(s[len('300BCG'):] if isinstance(s, str) else s) for s in samples])


def hash_sample_ids(samples):
    samples = make_integer_samples_ids(samples)
    return hash('_'.join(map(str, np.unique(samples))))


def _L2(estimator, X=None, y=None, sample_weight=None):
    return np.sum(np.square(estimator.named_steps.estimator.coef_))


def overpenalized_grid_cv(best_score_fraction, estimator, X, y, groups, cv, param_grid, scorer, scale=True,
                          feature_selection_score=None, n_features=None, n_pcs=None, grid_train_scores=False,
                          n_jobs=None, verbose=1):

    assert isinstance(estimator, LogisticRegression)
    assert not isinstance(estimator, RidgeClassifierCV)
    assert len(cv) < len(X)
    grid_search = grid_cv(estimator, X, y, groups, cv, param_grid, {'score': scorer, 'L2': _L2}, scale=scale,
                          feature_selection_score=feature_selection_score, n_features=n_features, n_pcs=n_pcs,
                          grid_train_scores=grid_train_scores, refit=False, n_jobs=n_jobs, verbose=verbose)

    cv_df = pd.DataFrame.from_dict(grid_search.cv_results_)
    max_param_idx = cv_df['mean_test_score'].idxmax()
    best_param_idx = cv_df.where(cv_df['mean_test_score'] > cv_df['mean_test_score'].max() * best_score_fraction)['mean_test_L2'].idxmin()
    print('Max. score: {:.3f} with {} L2 {}'.format(cv_df.loc[max_param_idx, 'mean_test_score'], cv_df.loc[max_param_idx, 'params'], cv_df.loc[max_param_idx, 'mean_test_L2']))
    print('Best LR_L2: {:.3f} with {} L2 {}'.format(cv_df.loc[best_param_idx, 'mean_test_score'], cv_df.loc[best_param_idx, 'params'], cv_df.loc[best_param_idx, 'mean_test_L2']))

    # You can test this by setting:
    # best_param_idx = max_param_idx

    grid_search.refit = True
    grid_search.best_index_ = best_param_idx
    grid_search.best_params_ = cv_df.loc[best_param_idx, 'params']
    grid_search.estimator.set_params(**grid_search.best_params_)
    grid_search.estimator.fit(X, y)
    grid_search.best_estimator_ = grid_search.estimator
    grid_search.best_score_ = cv_df.loc[best_param_idx, 'mean_test_score']

    return grid_search


def _make_score_func(scorer, y_true, y_pred, binarized=None, sample_weight=None):
    if isinstance(scorer, BinarizedRegressionScorer):
        y_true, y_pred = _binarize(y=y_true, y_pred=y_pred, binarized=binarized)
    if sample_weight:
        return scorer._sign * scorer._score_func(y_true, y_pred, sample_weight=sample_weight, **scorer._kwargs)
    else:
        return scorer._sign * scorer._score_func(y_true, y_pred, **scorer._kwargs)


def get_scorer_func(scorer, with_sample_weight=False, binarized=None):
    if with_sample_weight:
        scorer_func = lambda y_true, y_pred, sample_weight: _make_score_func(
            scorer=scorer, y_true=y_true, y_pred=y_pred, binarized=binarized, sample_weight=sample_weight)
    else:
        scorer_func = lambda y_true, y_pred: _make_score_func(
            scorer=scorer, y_true=y_true, y_pred=y_pred, binarized=binarized)
    return scorer_func


def frozen_nested_cv(estimator, X, y, y_real, X_extra, y_extra, groups, samples, target, test_kfold, cv_kfold, param_grid, scoring, scale=True,
                     feature_selection_name=None, feature_selection_mask=None, feature_selection_score=None,
                     n_features=None, add_extra_samples_for_fs=False, use_y_real_for_fs=False, n_pcs=None,
                     binary_mapper=None, grid_train_scores=False, alt_scoring=None, n_jobs=None, random_state=None, verbose=1, out=sys.stdout):

    grids, y_pred, y_true, samples_pred, feature_importances, selected_features, train_scores, test_scores, cv_scores =\
        [], [], [], [], [], [], [], [], []

    # to ensure that internal CV will have variability if setup with dict(shuffle=True, random_state=None)
    # especially in the case with test CV is leave-one-out
    if random_state is not None:
        print('Setting random seed to {}'.format(random_state), file=out)
        np.random.seed(random_state)

    print('', file=out)
    for k, (train, test) in enumerate(test_kfold.split(
            X, binary_mapper(y) if binary_mapper is not None else y, groups)):
        X_train, y_train = X[train], y[train]
        X_test, y_test, samples_test = X[test], y[test], samples[test]
        assert groups is None or len(set(groups[train]).intersection(groups[test])) == 0
        groups_train = groups[train] if groups is not None else None

        scorer = get_scorer(scoring)
        if binary_mapper is not None:
            scorer = binarize_labels(scorer, binary_mapper)

        if n_features is not None:
            if n_features >= X.shape[1]:
                add_masked_features = False
                feature_idx = np.arange(X.shape[1])
            else:
                add_masked_features = True
                if add_extra_samples_for_fs:
                    fs_X = np.concatenate([X_train[:, feature_selection_mask],
                                           X_extra[:, feature_selection_mask]], axis=0)
                    fs_Y = np.concatenate([y_real[train] if use_y_real_for_fs else y_train,
                                           y_extra], axis=0)
                else:
                    fs_X = X_train[:, feature_selection_mask]
                    fs_Y = y_real[train] if use_y_real_for_fs else y_train
                feature_scores = feature_selection_score(fs_X, fs_Y, k)
                assert feature_scores.shape[0] == X.shape[1]

                if 'blind_DE' in feature_selection_name:
                    # just a sanity check
                    from misc import de_fn
                    _is_LOO = isinstance(test_kfold, LeaveOneOut) or isinstance(test_kfold, LeaveOneGroupOut)
                    fn = de_fn(celltype='PBMC', model='{}.{}{}.seed{}_k{}'.format(
                        'LOO_samples' if _is_LOO else 'binarized_samples' if 'BIN_blind_DE' in feature_selection_name else 'test_samples',
                        target, '_score' if 'SCORE_blind_DE' in feature_selection_name else '',
                        100 if _is_LOO else random_state, k + 1), data='design')
                    de_samples = pd.read_csv(fn, index_col=0).index
                    de_donors = de_samples.str.replace('_V[1-3]_PBMC', '')
                    assert not np.isin(samples_test, de_samples.values).any()
                    assert not np.isin(samples_test, de_donors.values).any()

                if n_features >= 1:
                    feature_idx = np.argsort(feature_scores)[-n_features:]
                else:
                    feature_idx = np.where(feature_scores != 0)[0]
                X_train = np.concatenate([X_train[:, feature_idx], X_train[:, ~feature_selection_mask]], axis=1)
                X_test = np.concatenate([X_test[:, feature_idx], X_test[:, ~feature_selection_mask]], axis=1)

        _grid_cv = [(cv_train, cv_test) for cv_train, cv_test in cv_kfold.split(
            X_train, binary_mapper(y_train) if binary_mapper is not None else y_train, groups_train)]
        _grid_search_args = dict(estimator=estimator, X=X_train, y=y_train,
                                    groups=groups_train, cv=_grid_cv,
                                    param_grid=param_grid, scorer=scorer, scale=scale,
                                    feature_selection_score=None,
                                    n_features=None, n_pcs=n_pcs, grid_train_scores=grid_train_scores,
                                    n_jobs=n_jobs, verbose=verbose)
        if feature_selection_name and feature_selection_name.startswith('OP_'):
            grid_search = overpenalized_grid_cv(best_score_fraction=0.90, **_grid_search_args)
        else:
            grid_search = grid_cv(**_grid_search_args)

        # Test with the best estimator
        _model = grid_search.best_estimator_
        _estimator = _model.named_steps.estimator if isinstance(_model, Pipeline) else _model

        if hasattr(_model, 'decision_function'):
            _y_pred = _model.decision_function(X_test)
            _y_train_pred = _model.decision_function(X_train)
        elif hasattr(_model, 'predict_proba'):
            _y_pred = _model.predict_proba(X_test)[:, _model.classes_ == 1].flatten()
            _y_train_pred = _model.predict_proba(X_train)[:, _model.classes_ == 1].flatten()
        else:
            _y_pred = _model.predict(X_test)
            _y_train_pred = _model.predict(X_train)

        if isinstance(test_kfold, LeaveOneOut) or isinstance(test_kfold, LeaveOneGroupOut):
            test_scoring = 'ACC' if is_classifier(_estimator) else 'MSE'
            _loo_scoring = accuracy_score if is_classifier(_estimator) else mean_squared_error
            _test_score = _loo_scoring(binary_mapper(y_test) if binary_mapper is not None else y_test, _model.predict(X_test))
        else:
            test_scoring = scoring
            try:
                _test_score = get_scorer_func(scorer)(binary_mapper(y_test) if binary_mapper is not None else y_test, _y_pred)
            except ValueError as e:
                _test_score = np.nan
                print('WARNING:', str(e), file=out)

        _train_score = get_scorer_func(scorer)(binary_mapper(y_train) if binary_mapper is not None else y_train, _y_train_pred)

        grids.append(grid_search)
        y_true.append(binary_mapper(y_test) if binary_mapper is not None else y_test)
        y_pred.append(_y_pred)
        samples_pred.append(samples_test)
        feature_importances.append(_estimator.feature_importances_.reshape(-1,) if hasattr(_estimator, 'feature_importances_') \
                                       else (_estimator.coef_.reshape(-1,) if hasattr(_estimator, 'coef_') \
                                             else None))
        selected_features.append(
            (np.concatenate([feature_idx, np.where(~feature_selection_mask)[0]], axis=0) if add_masked_features else feature_idx)\
                if n_features is not None else None
        )
        test_scores.append(_test_score)
        train_scores.append(_train_score)
        cv_scores.append(grid_search.best_score_)

        if verbose > 0:
            if grid_train_scores:
                best_params_mask = np.array(grid_search.cv_results_['params']) == grid_search.best_params_
                assert best_params_mask.sum() == 1
                train_cv_score = grid_search.cv_results_['mean_train_score'][best_params_mask][0]
            print(
                'CV_{k} {scoring} {cv_score:.3f}{train_cv_score} {params}{coefs} TEST_{k} {test_scoring} {test_score:.3f} (train {scoring} {train_score:.3f}){n_iter}'.format(
                    k=k, scoring=scoring, test_scoring=test_scoring,
                    cv_score=grid_search.best_score_,
                    train_cv_score=' (train {} {:.3f})'.format(scoring, train_cv_score) if grid_train_scores else '',
                    params={
                        p[len(ESTIMATOR_PREFIX):] if p.startswith(ESTIMATOR_PREFIX) else p: grid_search.best_params_[p]
                        for p in grid_search.best_params_},
                    coefs=' coefs {}/{}'.format((_estimator.coef_ != 0).sum(),
                                                _estimator.coef_.flatten().shape[0]) if hasattr(_estimator,
                                                                                                'coef_') else '',
                    test_score=_test_score, train_score=_train_score,
                    n_iter=' n_iter {}{}'.format(_estimator.n_iter_,
                                                 '/{}'.format(_estimator.max_iter) if hasattr(_estimator, 'max_iter') else '') \
                        if hasattr(_estimator, 'n_iter_') else '',
                ),
                file=out
            )
            if grid_train_scores:
                print('mean_train_score',
                      ' '.join(list(map(lambda x: '{:<6.3f}'.format(x), grid_search.cv_results_['mean_train_score']))), file=out)
                print('mean_test_score',
                      ' '.join(list(map(lambda x: '{:<6.3f}'.format(x), grid_search.cv_results_['mean_test_score']))), file=out)
            out.flush()

    if verbose > 0:
        print('\nCV_mean {scoring} {cv_score:.3f} TEST_mean {test_scoring} {test_score:.3f} (train_mean {scoring} {train_score:.3f})'.format(
            scoring=scoring, test_scoring=test_scoring,
            cv_score=np.mean(cv_scores), test_score=np.mean(test_scores), train_score=np.mean(train_scores)), file=out)
        print(
            'TEST_union {scoring} {test_score}'.format(
                scoring=scoring, test_score=get_scorer_func(scorer)(
                    np.asarray(binary_mapper([v for _y in y_true for v in _y]) if binary_mapper is not None else [v for _y in y_true for v in _y]),
                    np.asarray([v for _y in y_pred for v in _y]))),
                file=out)
        if alt_scoring is not None:
            print(
                'TEST_union {scoring} {test_score}'.format(
                    scoring=alt_scoring, test_score=get_scorer_func(get_scorer(alt_scoring))(
                        np.asarray(binary_mapper([v for _y in y_true for v in _y]) if binary_mapper is not None else [v for _y in y_true for v in _y]),
                        np.asarray([v for _y in y_pred for v in _y]))),
                    file=out)
        out.flush()

    return grids, np.asarray(y_true), np.asarray(y_pred), np.asarray(samples_pred),\
           np.asarray(feature_importances), np.asarray(selected_features),\
           np.asarray(train_scores), np.asarray(test_scores)


def nested_cv_with_blind_DE(estimator, X, y, groups, samples, peaks, target, test_kfold, cv_kfold, param_grid, scoring, scale=True,
              feature_selection_score=None, n_features=None, random_features=False, n_pcs=None, binary_mapper=None, grid_train_scores=False,
              alt_scoring=None, n_jobs=None, random_state=None, verbose=1, out=sys.stdout):

    assert n_features is None or n_features < X.shape[1]

    # to ensure that internal CV will have variability if setup with dict(shuffle=True, random_state=None)
    # especially in the case with test CV is leave-one-out
    if random_state is not None:
        print('Setting random seed to {}'.format(random_state), file=out)
        np.random.seed(random_state)

    # This will make the hashing work
    X = np.concatenate([make_integer_samples_ids(samples).reshape((-1, 1)), X], axis=1)

    grids, y_pred, y_true, samples_pred, feature_importances, selected_features, train_scores, test_scores, cv_scores =\
        [], [], [], [], [], [], [], [], []

    for k_test, (train, test) in enumerate(test_kfold.split(
            X, binary_mapper(y) if binary_mapper is not None else y, groups)):
        X_train, y_train = X[train], y[train]
        X_test, y_test, samples_test = X[test], y[test], samples[test]
        assert groups is None or len(set(groups[train]).intersection(groups[test])) == 0
        groups_train = groups[train] if groups is not None else None

        scorer = get_scorer(scoring)
        if binary_mapper is not None:
            scorer = binarize_labels(scorer, binary_mapper)

        # just a sanity check
        from misc import de_fn, ML_RESULTS_ROOT
        fn = de_fn(celltype='PBMC', model='test_samples_nested_10_10.{}.seed{}_test{}'.format(
            target, random_state, k_test + 1), data='{data}', results_dir=ML_RESULTS_ROOT)
        de_samples = pd.read_csv(fn.format(data='design'), index_col=0).index
        de_donors = de_samples.str.replace('_V[1-3]_PBMC', '')
        assert not np.isin(samples[test], de_samples.values).any()
        assert not np.isin(samples[test], de_donors.values).any()
        de_peaks = pd.read_csv(fn.format(data='results_p5'), usecols=[0], squeeze=True).index.values
        assert np.array_equal(peaks, de_peaks)

        _grid_cv = [(cv_train, cv_test) for cv_train, cv_test in cv_kfold.split(
            X_train, binary_mapper(y_train) if binary_mapper is not None else y_train, groups_train)]

        # just a sanity check
        for k_cv, (cv_train, cv_test) in enumerate(_grid_cv):
            fn = de_fn(celltype='PBMC', model='test_samples_nested_10_10.{}.seed{}_test{}_cv{}'.format(
                target, random_state, k_test + 1, k_cv + 1), data='{data}', results_dir=ML_RESULTS_ROOT)
            de_samples = pd.read_csv(fn.format(data='design'), index_col=0).index
            de_donors = de_samples.str.replace('_V[1-3]_PBMC', '')
            assert not np.isin(samples[train][cv_test], de_samples.values).any()
            assert not np.isin(samples[train][cv_test], de_donors.values).any()
            de_peaks = pd.read_csv(fn.format(data='results_p5'), usecols=[0], squeeze=True).index.values
            assert np.array_equal(peaks, de_peaks)

        grid_search = grid_cv(estimator=estimator, X=X_train, y=y_train, groups=groups_train, cv=_grid_cv,
                              param_grid=param_grid, scorer=scorer, scale=scale,
                              feature_selection_score=feature_selection_score,
                              n_features=n_features, random_features=random_features,
                              n_pcs=n_pcs, grid_train_scores=grid_train_scores,
                              n_jobs=n_jobs, verbose=verbose)

        # Test with the best estimator
        _model = grid_search.best_estimator_
        _estimator = _model.named_steps.estimator

        if hasattr(_model, 'decision_function'):
            _y_pred = _model.decision_function(X_test)
            _y_train_pred = _model.decision_function(X_train)
        elif hasattr(_model, 'predict_proba'):
            _y_pred = _model.predict_proba(X_test)[:, _model.classes_ == 1].flatten()
            _y_train_pred = _model.predict_proba(X_train)[:, _model.classes_ == 1].flatten()
        else:
            _y_pred = _model.predict(X_test)
            _y_train_pred = _model.predict(X_train)
        _test_score = get_scorer_func(scorer)(binary_mapper(y_test) if binary_mapper is not None else y_test, _y_pred)
        _train_score = get_scorer_func(scorer)(binary_mapper(y_train) if binary_mapper is not None else y_train, _y_train_pred)

        grids.append(grid_search)
        y_true.append(binary_mapper(y_test) if binary_mapper is not None else y_test)
        y_pred.append(_y_pred)
        samples_pred.append(samples_test)
        feature_importances.append(_estimator.feature_importances_.reshape(-1,) if hasattr(_estimator, 'feature_importances_') \
                                       else (_estimator.coef_.reshape(-1,) if hasattr(_estimator, 'coef_') \
                                             else None))
        if n_features is not None:
            _selected = _model.named_steps.feature_selector.get_support()
            assert _selected.dtype == bool
            _selected = _selected[1:]
        selected_features.append(_selected if n_features is not None else None)
        test_scores.append(_test_score)
        train_scores.append(_train_score)
        cv_scores.append(grid_search.best_score_)

        if verbose > 0:
            if grid_train_scores:
                best_params_mask = np.array(grid_search.cv_results_['params']) == grid_search.best_params_
                assert best_params_mask.sum() == 1
                train_cv_score = grid_search.cv_results_['mean_train_score'][best_params_mask][0]
            print(
                'CV_{k} {scoring} {cv_score:.3f}{train_cv_score} {params}{coefs} TEST_{k} {scoring} {test_score:.3f} (train {scoring} {train_score:.3f}){n_iter}'.format(
                    k=k_test, scoring=scoring, cv_score=grid_search.best_score_,
                    train_cv_score=' (train {} {:.3f})'.format(scoring, train_cv_score) if grid_train_scores else '',
                    params={
                        p[len(ESTIMATOR_PREFIX):] if p.startswith(ESTIMATOR_PREFIX) else p: grid_search.best_params_[p]
                        for p in grid_search.best_params_},
                    coefs=' coefs {}/{}'.format((_estimator.coef_ != 0).sum(),
                                                _estimator.coef_.flatten().shape[0]) if hasattr(_estimator,
                                                                                                'coef_') else '',
                    test_score=_test_score, train_score=_train_score,
                    n_iter=' n_iter {}{}'.format(_estimator.n_iter_,
                                                 '/{}'.format(_estimator.max_iter) if hasattr(_estimator, 'max_iter') else '') \
                        if hasattr(_estimator, 'n_iter_') else '',
                ),
                file=out
            )
            if grid_train_scores:
                print('mean_train_score',
                      ' '.join(list(map(lambda x: '{:<6.3f}'.format(x), grid_search.cv_results_['mean_train_score']))), file=out)
                print('mean_test_score',
                      ' '.join(list(map(lambda x: '{:<6.3f}'.format(x), grid_search.cv_results_['mean_test_score']))), file=out)

    if verbose > 0:
        print('CV_mean {scoring} {cv_score:.3f} TEST_mean {scoring} {test_score:.3f} (train_mean {scoring} {train_score:.3f})'.format(
            scoring=scoring, cv_score=np.mean(cv_scores), test_score=np.mean(test_scores), train_score=np.mean(train_scores)), file=out)
        print(
            'TEST_union {scoring} {test_score}'.format(
                scoring=scoring, test_score=get_scorer_func(scorer)(
                    np.asarray(binary_mapper([v for _y in y_true for v in _y]) if binary_mapper is not None else [v for _y in y_true for v in _y]),
                    np.asarray([v for _y in y_pred for v in _y]))),
                file=out)
        if alt_scoring is not None:
            print(
                'TEST_union {scoring} {test_score}'.format(
                    scoring=alt_scoring, test_score=get_scorer_func(get_scorer(alt_scoring))(
                        np.asarray(binary_mapper([v for _y in y_true for v in _y]) if binary_mapper is not None else [v for _y in y_true for v in _y]),
                        np.asarray([v for _y in y_pred for v in _y]))),
                    file=out)

    return grids, np.asarray(y_true), np.asarray(y_pred), np.asarray(samples_pred),\
           np.asarray(feature_importances), np.asarray(selected_features),\
           np.asarray(train_scores), np.asarray(test_scores)


def _binarize(y, y_pred, binarized):
    mask = (y <= binarized[0]) | (y >= binarized[1])
    binary_y = np.copy(y[mask])
    binary_y_pred = np.copy(y_pred[mask])
    binary_y[y[mask] <= binarized[0]] = 0
    binary_y[y[mask] >= binarized[1]] = 1
    assert set(binary_y).issubset([0, 1])
    assert len(binary_y) == len(binary_y_pred)
    return binary_y, binary_y_pred


class SortOfStratifiedKFold(StratifiedKFold):
    def __init__(self, n_splits=5, *, shuffle=False, thresholds=None, random_state=None):
        super().__init__(n_splits=n_splits, shuffle=shuffle, random_state=random_state)
        self._thresholds = sorted(thresholds)
        print(self._thresholds)

    def split(self, X, y, groups=None):
        binary_y = np.copy(y)
        previous_thr = -np.inf
        for i, thr in enumerate(self._thresholds):
            binary_y[(y > previous_thr) & (y <= thr)] = i
            previous_thr = thr
        binary_y[y >= self._thresholds[-1]] = len(self._thresholds)
        return super().split(X, binary_y, groups)


class BinarizedRegressionScorer(_BaseScorer):
    def __init__(self, score_func, sign, binarized):
        super(BinarizedRegressionScorer, self).__init__(score_func=score_func, sign=sign, kwargs=dict())
        self._binarized = binarized

    def _score(self, method_caller, estimator, X, y, sample_weight=None):
        y_pred = method_caller(estimator, "predict", X)
        binary_y, binary_y_pred = _binarize(y=y, y_pred=y_pred, binarized=self._binarized)
        if sample_weight is not None:
            return self._sign * self._score_func(binary_y, binary_y_pred, sample_weight=sample_weight, **self._kwargs)
        else:
            return self._sign * self._score_func(binary_y, binary_y_pred, **self._kwargs)


def roc_auc_score25(y_true, y_score):
    return roc_auc_score(y_true=y_true, y_score=y_score, max_fpr=0.25)


def nested_cv(estimator, X, y, groups, samples, target, test_kfold, cv_kfold, param_grid, scoring,
              permute_labels=False, subsample=False, scale=True,
              feature_selection_score=None, n_features=None, n_pcs=None, pca_first=False, binary_mapper=None, binarized=None,
              merged_folds_scoring=False, grid_train_scores=False, alt_scoring=None, strictly_use_proba=False,
              n_jobs=None, random_state=None, verbose=1, out=sys.stdout):

    print('subsample', subsample)
    # to ensure that internal CV will have variability if setup with dict(shuffle=True, random_state=None)
    # especially in the case with test CV is leave-one-out
    if random_state is not None:
        print('Setting random seed to {}'.format(random_state), file=out)
        np.random.seed(random_state)

    if subsample:
        subsample_rs = np.random.RandomState(random_state)
        subsample_idx = np.arange(len(y))
        if is_classifier(estimator):
            assert set(y) in [set([0, 1]), set([-1, 1])]
            idx1 = subsample_rs.choice(subsample_idx[y == 1], size=int(np.ceil(len(subsample_idx[y == 1]) * subsample)), replace=False)
            idx2 = subsample_rs.choice(subsample_idx[y != 1], size=int(np.ceil(len(subsample_idx[y != 1]) * subsample)), replace=False)
            subsample_idx = np.concatenate([idx1, idx2])
        else:
            subsample_idx = subsample_rs.choice(subsample_idx, size=int(np.ceil(len(subsample_idx) * subsample)), replace=False)
        assert len(set(subsample_idx)) == len(subsample_idx) and len(subsample_idx) <= len(y)
        print('Subsampling', y.shape, '-->', end=' ')
        print(subsample_idx.shape)
    else:
        subsample_idx = None

    if permute_labels:
        permutation_rs = np.random.RandomState(random_state)

    grids, y_pred, y_true, samples_pred, feature_importances, selected_features, train_scores, test_scores, cv_scores =\
        [], [], [], [], [], [], [], [], []
    print(test_kfold)
    print(cv_kfold)
    for k_test, (train, test) in enumerate(test_kfold.split(X, binary_mapper(y) if binary_mapper is not None else y, groups)):

        if subsample:
            assert np.isin(train, np.arange(len(y))).all()
            train = train[np.isin(train, subsample_idx)]

        if permute_labels:
            assert np.isin(train, np.arange(len(y))).all()
            permutation_idx = permutation_rs.choice(len(train), size=len(train), replace=False)
            permuted_train = train[permutation_idx]

        X_train = X[train]
        y_train = y[train] if not permute_labels else y[permuted_train]
        X_test, y_test, samples_test = X[test], y[test], samples[test]
        assert groups is None or len(set(groups[train]).intersection(groups[test])) == 0
        groups_train = groups[train] if groups is not None else None
        if scoring in ['quantile_roc_auc', 'quantile_average_precision']:
            assert binary_mapper is None
            assert binarized is not None
            if scoring == 'quantile_roc_auc':
                score_func = roc_auc_score
            elif scoring == 'quantile_average_precision':
                score_func = average_precision_score
            else:
                raise ValueError
            scorer = BinarizedRegressionScorer(score_func, sign=1, binarized=binarized)
        elif scoring == 'partial_roc_auc':
            scorer = make_scorer(score_func=roc_auc_score25, greater_is_better=True, needs_threshold=True)
        else:
            scorer = get_scorer(scoring)
            if binary_mapper is not None:
                scorer = binarize_labels(scorer, binary_mapper)

        _grid_cv = [(cv_train, cv_test) for cv_train, cv_test in cv_kfold.split(
            X_train, binary_mapper(y_train) if binary_mapper is not None else y_train, groups_train)]

        if scale and isinstance(estimator, (BaggingClassifier, RidgeClassifierCV)):
            scaler = StandardScaler()
            X_train = scaler.fit_transform(X_train)
            X_test = scaler.transform(X_test)

        grid_search = grid_cv(estimator=estimator, X=X_train, y=y_train, groups=groups_train, cv=_grid_cv,
                              param_grid=param_grid, scorer=scorer,
                              scale=scale and not isinstance(estimator, (BaggingClassifier, RidgeClassifierCV)),
                              feature_selection_score=feature_selection_score,
                              n_features=n_features, n_pcs=n_pcs, pca_first=pca_first,
                              merged_folds_scoring=merged_folds_scoring, grid_train_scores=grid_train_scores,
                              n_jobs=n_jobs, verbose=verbose)

        # Test with the best estimator
        _model = grid_search.best_estimator_
        _estimator = _model.named_steps.estimator if isinstance(_model, Pipeline) else _model

        if not strictly_use_proba and hasattr(_model, 'decision_function'):
            _y_pred = _model.decision_function(X_test)
            _y_train_pred = _model.decision_function(X_train)
        elif hasattr(_model, 'predict_proba'):
            _y_pred = _model.predict_proba(X_test)[:, _model.classes_ == 1].flatten()
            _y_train_pred = _model.predict_proba(X_train)[:, _model.classes_ == 1].flatten()
        else:
            _y_pred = _model.predict(X_test)
            _y_train_pred = _model.predict(X_train)

        if len(y_test) < 3 and not isinstance(test_kfold, (LeaveOneOut, LeaveOneGroupOut)):
            print('WARNING: Less than three test samples', file=out)
        _test_score = np.nan if len(y_test) < 3 or (is_classifier(estimator) and len(set(y_test)) < 2) \
            else get_scorer_func(scorer, binarized=binarized)(binary_mapper(y_test) if binary_mapper is not None else y_test, _y_pred)
        _train_score = get_scorer_func(scorer, binarized=binarized)(binary_mapper(y_train) if binary_mapper is not None else y_train, _y_train_pred)

        grids.append(grid_search)
        y_true.append(binary_mapper(y_test) if binary_mapper is not None else y_test)
        y_pred.append(_y_pred)
        samples_pred.append(samples_test)
        feature_importances.append(_estimator.feature_importances_.reshape(-1,) if hasattr(_estimator, 'feature_importances_') \
                                       else (_estimator.coef_.reshape(-1,) if hasattr(_estimator, 'coef_') \
                                             else None))
        selected_features.append(_model.named_steps.feature_selector.get_support() if n_features is not None else None)
        test_scores.append(_test_score)
        train_scores.append(_train_score)
        cv_scores.append(grid_search.best_score_)

        if verbose > 0:
            if grid_train_scores:
                best_params_mask = np.array(grid_search.cv_results_['params']) == grid_search.best_params_
                assert best_params_mask.sum() == 1
                train_cv_score = grid_search.cv_results_['mean_train_score'][best_params_mask][0]
            print(
                'CV_{k} {scoring} {cv_score:.3f}{train_cv_score} {params}{coefs} TEST_{k} {scoring} {test_score:.3f} (train {scoring} {train_score:.3f}){n_iter}'.format(
                    k=k_test, scoring=scoring, cv_score=grid_search.best_score_,
                    train_cv_score=' (train {} {:.3f})'.format(scoring, train_cv_score) if grid_train_scores else '',
                    params={
                        p[len(ESTIMATOR_PREFIX):] if p.startswith(ESTIMATOR_PREFIX) else p: grid_search.best_params_[p]
                        for p in grid_search.best_params_},
                    coefs=' coefs {}/{}'.format((_estimator.coef_ != 0).sum(),
                                                _estimator.coef_.flatten().shape[0]) if hasattr(_estimator,
                                                                                                'coef_') else '',
                    test_score=_test_score, train_score=_train_score,
                    n_iter=' n_iter {}{}'.format(_estimator.n_iter_,
                                                 '/{}'.format(_estimator.max_iter) if hasattr(_estimator, 'max_iter') else '') \
                        if hasattr(_estimator, 'n_iter_') else '',
                ),
                file=out
            )
            if grid_train_scores:
                print('mean_train_score',
                      ' '.join(list(map(lambda x: '{:<6.3f}'.format(x), grid_search.cv_results_['mean_train_score']))), file=out)
                print('mean_test_score',
                      ' '.join(list(map(lambda x: '{:<6.3f}'.format(x), grid_search.cv_results_['mean_test_score']))), file=out)

    if verbose > 0:
        test_scores = np.asarray(test_scores)
        _result_msg = 'CV_mean {scoring} {cv_score:.3f} TEST_mean {scoring} {test_score:.3f} (train_mean {scoring} {train_score:.3f})'.format(
            scoring=scoring, cv_score=np.mean(cv_scores), test_score=np.mean(test_scores[~np.isnan(test_scores)]), train_score=np.mean(train_scores))
        print(_result_msg)
        print(_result_msg, file=out)
        print(
            'TEST_union {scoring} {test_score}'.format(
                scoring=scoring, test_score=get_scorer_func(scorer, binarized=binarized)(
                    np.asarray(binary_mapper([v for _y in y_true for v in _y]) if binary_mapper is not None else [v for _y in y_true for v in _y]),
                    np.asarray([v for _y in y_pred for v in _y]))),
                file=out)
        if alt_scoring is not None:
            print(
                'TEST_union {scoring} {test_score}'.format(
                    scoring=alt_scoring,
                    test_score=get_scorer_func(get_scorer(alt_scoring) if alt_scoring not in ['quantile_roc_auc', 'quantile_average_precision', 'partial_roc_auc'] else \
                                                   make_scorer(score_func=roc_auc_score25, greater_is_better=True, needs_threshold=True) if alt_scoring == 'partial_roc_auc' else \
                                                   BinarizedRegressionScorer(average_precision_score if alt_scoring == 'quantile_average_precision' else roc_auc_score, sign=1, binarized=binarized),
                                               binarized=binarized)(
                        np.asarray(binary_mapper([v for _y in y_true for v in _y]) if binary_mapper is not None else [v for _y in y_true for v in _y]),
                        np.asarray([v for _y in y_pred for v in _y])
                    )
                ), file=out)

    return grids, np.asarray(y_true), np.asarray(y_pred), np.asarray(samples_pred),\
           np.asarray(feature_importances), np.asarray(selected_features),\
           np.asarray(train_scores), np.asarray(test_scores)


def df_append(df, name=None, **kwargs):
    s = pd.Series(list(kwargs.values()), index=list(kwargs.keys()), name=name if name is not None else len(df))
    return df.append(s, sort=False)


def get_majority_class(y, positive_class=1, negative_class=0):
    assert isinstance(y, np.ndarray)
    assert list(np.unique(y[~np.isnan(y)])) == [negative_class, positive_class]
    majority_class = positive_class if (y == positive_class).sum() >= (y == negative_class).sum() else negative_class
    return majority_class


def run_prerank_gseapy(scores, gene_sets, min_size=0, max_size=20000, permutation_num=1000, weighted_score_type=1,
                       ascending=False, processes=8, seed=0, gene_sets_name=None, description=None,
                       outdir=None):
    import gseapy
    results = gseapy.prerank(rnk=scores, gene_sets=gene_sets, permutation_num=permutation_num,
                          weighted_score_type=weighted_score_type, outdir=outdir,
                          min_size=min_size, max_size=max_size, ascending=ascending, processes=processes,
                          seed=seed, no_plot=True)
    df = results.res2d.reset_index()
    if df.shape[0] > 0:
        df = df.drop(['geneset_size', 'genes'], axis=1)
        df = df.rename({'Term': 'description', 'es': 'ES', 'nes': 'NES', 'pval': 'p_value',
                        'fdr': 'adjusted_p_value', 'matched_size': 'gene_set_size'}, axis=1)
        df = df.reindex(df['NES'].abs().sort_values(ascending=False).index).reset_index(drop=True).rename_axis(
            'rank').reset_index(drop=False)
        df['rank'] = df['rank'] + 1
        if gene_sets_name is not None:
            df['gene_set_library'] = gene_sets_name
        if description is not None:
            df['comparison'] = description
        df['ledge_genes'] = df['ledge_genes'].str.split(';').map(str)

    return df, results


def read_gene_sets(fn, description_NA='--', as_arrays=True, assert_upper=True, make_upper=False):
    '''
    Reads in a gene set file in a GMT format.

    :param fn: gene set filename downloaded from Enrichr website
    :param description_NA: I assume description is missing or equal to "description_NA" symbol for safety
    :param as_arrays: whether the second item of the tuple should be an array or a list
    :param assert_upper: fail if there is a gene with lower-case characters
    :return: iterable of tuples (term, genes)
    '''

    with open_file(fn) as f:
        for line in f:
            line_split = line.strip().split('\t')
            term, description, genes = line_split[0].strip(), line_split[1].strip(), line_split[2:]
            assert len(term) != 0
            assert description_NA is None or len(description) == 0 or description == description_NA
            assert all([' ' not in g and ',' not in g for g in genes]), 'Illegal character " " or "," in a gene name'
            assert not assert_upper or all([g == g.upper() for g in genes]), 'Lower-case character in a gene name'
            if make_upper:
                genes = list(map(str.upper, genes))
            yield term, np.asarray(genes) if as_arrays else genes


def gene_set_library(fn, description_NA='--', as_arrays=True, assert_upper=True, make_upper=False):
    return {term: genes for term, genes in read_gene_sets(fn=fn, description_NA=description_NA, as_arrays=as_arrays,
                                                          assert_upper=assert_upper, make_upper=make_upper)}


def gene_set_enrichment_test(top_genes, gene_set_libraries, background, padj_method='fdr_bh', fisher_hypergeom_tol=1e-4,
                    check_hypergeom=False, check_background=False, adjust_gene_sets_to_background=False, padj_union=False,
                            min_gs=None, max_gs=None):
    assert len(top_genes) == len(set(top_genes))
    assert is_iterable(background) or not adjust_gene_sets_to_background
    assert is_iterable(background) or not check_background

    if check_background:
        assert all([g in background for g in top_genes])
        
    background_size = len(background) if is_iterable(background) else background

    res_df = pd.DataFrame(columns=ENRICHR_COLS)
    for lib_name, library in gene_set_libraries:
        results = []
        for term in library:
            gs = np.asarray(library[term]).copy()
            if adjust_gene_sets_to_background:
                gs = gs[np.isin(gs, background)]
            if check_background:
                assert all([g in background for g in gs])

            if min_gs is not None and len(gs) < min_gs:
                continue
                
            if max_gs is not None and len(gs) > max_gs:
                continue
            
            # hits are sorted in the order of top_genes
            hits = [g for g in top_genes if g in gs]
            overlap = '{}/{}'.format(len(hits), len(gs))

            if len(hits) != 0:
                term_not_top = set(gs).difference(top_genes)
                top_not_term = set(top_genes).difference(gs)
                term_or_top = set(gs).union(top_genes)
                oddsratio, f_p = stats.fisher_exact(
                    [[len(hits), len(top_not_term)], [len(term_not_top), background_size - len(term_or_top)]],
                    alternative='greater'
                )
                results.append((lib_name, term, f_p, np.nan, np.nan, oddsratio, overlap, ';'.join(hits)))

                if check_hypergeom:
                    k = len(hits)
                    M = background_size
                    n = len(gs)
                    N = len(top_genes)
                    hg_p = stats.hypergeom.sf(k - 1, M, n, N)
                    assert f_p - hg_p < fisher_hypergeom_tol

        if len(results) != 0:
            df = pd.DataFrame(results, columns=ENRICHR_COLS)
            df['Adjusted P-value'] = multipletests(df['P-value'].values, method=padj_method)[1]
            res_df = pd.concat([res_df, df.sort_values('P-value')])
        else:
            print('WARNING: No overlap for any term:', lib_name)

    if padj_union:
        res_df['Adjusted P-value'] = multipletests(res_df['P-value'].values, method=padj_method)[1]

    return res_df.reset_index(drop=True)


def page_test(LFCs, gene_set_mask, min_genes=None, max_genes=None, two_sided=True, ddof=0, axis=1):
    '''
    Performs the PAGE test for one gene set but possible across many comparisons.

    :param LFCs: a dataframe of log2 fold changes; when axis=1, genes are as columns and comparisons as rows
    :param gene_set_mask: a set of genes or a boolean mask for the genes in the LFCs dataframe
    :param min_genes: do not test if less than min_genes would be tested for the given gene set
    :param max_genes: do not test if more than max_genes would be tested for the given gene set
    :param two_sided: two-sided T-test
    :param ddof: degrees of freedom
    :param axis: 1 if genes are in the columns and 0 otherwise
    :return: a tuple of arrays with a value for each comparison (z_score, p_value, mean, std,
    gene_set_mean, gene_set_sizes, gene_set_mask)
    '''

    if axis == 1:
        # genes need to be along the axis 0
        LFCs = LFCs.T

    if gene_set_mask.dtype != bool:
        gene_set_mask = LFCs.index.isin(gene_set_mask)

    if (min_genes is None or gene_set_mask.sum() >= min_genes) \
            and (max_genes is None or gene_set_mask.sum() <= max_genes):
        gene_set_mask = np.broadcast_to(gene_set_mask, (LFCs.shape[1], gene_set_mask.shape[0])).T
        gene_set_mask = gene_set_mask & ~np.isnan(LFCs)
        gene_set_sizes = gene_set_mask.sum(axis=0)

        mean = LFCs.mean(axis=0)
        std = np.sqrt(LFCs.var(axis=0, ddof=ddof) / gene_set_sizes)
        gene_set_mean = LFCs[gene_set_mask].mean(axis=0)
        z = (gene_set_mean - mean) / std
        p = stats.norm.sf(np.abs(z))

        return (z, (p * 2) if two_sided else p, mean, std, gene_set_mean,
                gene_set_sizes, gene_set_mask.T if axis == 1 else gene_set_mask)
    else:
        return None, None, None, None, None, gene_set_mask.sum(), None


def page(LFCs, gene_set_libraries, gs_path=None, min_genes=None, max_genes=None, adjust_pvals=True, rank=True):
    '''
    Runs page tests across several gene set libraries and comparisons.

    :param LFCs: a dataframe of log2 fold changes; when axis=1, genes are as columns and comparisons as rows
    :param gene_set_libraries: names of the gene set libraries files
    :param gs_path: directory containing gene set libraries files
    :param min_genes: do not test if less than min_genes would be tested for the given gene set
    :param max_genes: do not test if more than max_genes would be tested for the given gene set
    :param adjust_pvals: perform Benjamini-Hochberg FDR correction
    :param rank: rank by absolute z-score
    :return: a long dataframe of all comparisons, gene set libraries, and terms
    '''

    if isinstance(gene_set_libraries, str):
        gene_set_libraries = [gene_set_libraries]

    page_df = pd.DataFrame()

    for gene_set_library in gene_set_libraries:
        lib_df = pd.DataFrame()
        for term, genes in read_gene_sets(os.path.join(gs_path, gene_set_library)):
            z_scores, p_values, mean, std, gene_set_mean, gene_set_sizes, gene_set_mask = \
                page_test(LFCs, genes, min_genes=min_genes, max_genes=max_genes)
            if z_scores is not None:
                term_df = pd.DataFrame(index=pd.MultiIndex.from_product(
                    [[gene_set_library], [term], z_scores.index],
                    names=['gene_set_library', 'description', 'comparison']
                ))
                term_df['z_score'] = z_scores.values
                term_df['p_value'] = p_values
                term_df['mean'] = mean.values
                term_df['std'] = std.values
                term_df['gene_set_mean'] = gene_set_mean.values
                term_df['delta_mean'] = (gene_set_mean - mean).values
                term_df['gene_set_size'] = gene_set_sizes.values
                assert constant(gene_set_mask.values, axis=0).all()

                # I was considering defining lead genes as all those that have LFCs larger than the mean
                # but I decided it against it in the end
                # lead_mask = (((LFCs.where(z_scores > 0).T >= np.clip(mean, a_min=0, a_max=None)) | (
                #             LFCs.where(z_scores < 0).T <= np.clip(mean, a_min=None, a_max=0))) & gene_set_mask.T).T

                all_genes, contrib_genes, lead_genes = [], [], []
                for c in z_scores.index:
                    gene_set_LFCs = LFCs.loc[c, gene_set_mask.loc[c]].sort_values(ascending=z_scores.loc[c] < 0)
                    all_genes.append(str(gene_set_LFCs.index.tolist()))
                    contrib_genes.append(str(gene_set_LFCs.loc[(gene_set_LFCs > 0) if z_scores.loc[c] > 0 else (gene_set_LFCs < 0)].index.tolist()))

                    # I was considering defining lead genes as all those that have LFCs larger than the mean
                    # but I decided it against it in the end
                    # lead_LFCs = LFCs.loc[c, lead_mask.loc[c]].sort_values(ascending=z_scores.loc[c] < 0)

                    # lead genes are those with within gene set absolute z-score more than 1
                    # Note: my personal definition
                    lead_LFCs = gene_set_LFCs
                    z_lead_LFCs = (lead_LFCs - lead_LFCs.mean()) / lead_LFCs.std()
                    lead_genes.append(str(lead_LFCs.loc[(z_lead_LFCs >= 1) if z_scores.loc[c] > 0 else (z_lead_LFCs <= -1)].index.tolist()))
                term_df['ledge_genes'] = lead_genes  # leading edge genes
                term_df['contrib_genes'] = contrib_genes  # gene set genes that were up/down if the pathway was up/down
                term_df['genes'] = all_genes  # all the genes from the gene set that were measured
                lib_df = pd.concat([lib_df, term_df], axis=0)

        if adjust_pvals:
            lib_df['adjusted_p_value'] = np.nan
            pval_mask = ~pd.isnull(lib_df['p_value'])
            for comparison in np.unique(lib_df.index.get_level_values('comparison')):
                cmp_mask = lib_df.index.get_level_values('comparison') == comparison
                lib_df.loc[cmp_mask & pval_mask, 'adjusted_p_value'] = \
                    multipletests(lib_df.loc[cmp_mask & pval_mask, 'p_value'].values, method='fdr_bh')[1]

        if rank:
            lib_df['rank'] = np.nan
            for comparison in np.unique(lib_df.index.get_level_values('comparison')):
                cmp_mask = lib_df.index.get_level_values('comparison') == comparison
                lib_df.loc[cmp_mask, 'rank'] = rankdata(-lib_df.loc[cmp_mask, 'z_score'].abs().values, method='ordinal')

        page_df = pd.concat([page_df, lib_df], axis=0)

    return page_df


def reduce_regions_to_genes(group):
    if (group > 0).sum() / len(group) >= 2 / 3:
        return group.max()
    elif (group < 0).sum() / len(group) >= 2 / 3:
        return group.min()
    else:
        return group.mean()


def kinship(X, normalize=True):
    if X.dtype != np.float64:
        X = np.array(X, dtype=np.float64)
    K = X.dot(X.T)
    if normalize:
        K /= K.diagonal().mean()
    return K


def compare_linear_models(ssr, df_resid, scale=None):
    assert len(ssr) == len(df_resid) & len(ssr) > 1
    assert sorted(df_resid)[::-1] == list(df_resid), list(df_resid)
    names = ['df_resid', 'ssr', 'df_diff', 'ss_diff', 'F', 'p']
    table = pd.DataFrame()
    table["ssr"] = ssr
    table["df_resid"] = df_resid
    if scale is None:
        if len(ssr) > 2:
            print('Warning: using scale of the last most model for all comparisons')
        scale = table["ssr"].iloc[-1] / table["df_resid"].iloc[-1]
    table["df_diff"] = -table["df_resid"].diff()
    table["ss_diff"] = -table["ssr"].diff()
    table["F"] = table["ss_diff"] / table["df_diff"] / scale
    mask = np.arange(len(ssr)) > 0
    table.loc[mask, "p"] = stats.f.sf(table["F"].loc[mask], table["df_diff"].loc[mask], table["df_resid"].loc[mask])
    return table


def adata_scatter_plot(adata, plot, color, ax=None, title=None, annotation=None,
                 figure_fn=None, set_rasterized=False, legend_title_fontsize=None, legend_kwargs=None, **kwargs):
    import scanpy as sc
    if plot == 'pca':
        plot_func = sc.pl.pca
    elif plot == 'tsne':
        plot_func = sc.pl.tsne
    elif plot == 'umap':
        plot_func = sc.pl.umap
    else:
        raise ValueError

    ax = plot_func(adata, color=color, ax=ax, show=False, **kwargs)
    if title is not None:
        ax.set_title(title)
    if legend_kwargs is not None:
        ax.legend(**legend_kwargs)
    if legend_title_fontsize is not None:
        ax.get_legend().get_title().set_fontsize(legend_title_fontsize)
    ax.set_rasterized(set_rasterized)
    sns.despine()

    if annotation is not None:
        annotate_scatter(adata.obs[annotation], adata.obsm['X_{}'.format(plot)][:, :2], ax)

    if figure_fn is not None:
        plt.savefig(figure_fn, bbox_inches='tight')

    return ax


# https://github.com/yandexdataschool/roc_comparison/blob/master/compare_auc_delong_xu.py


def _compute_midrank(x):
    """Computes midranks.
    Args:
       x - a 1D numpy array
    Returns:
       array of midranks
    """
    J = np.argsort(x)
    Z = x[J]
    N = len(x)
    T = np.zeros(N, dtype=np.float)
    i = 0
    while i < N:
        j = i
        while j < N and Z[j] == Z[i]:
            j += 1
        T[i:j] = 0.5*(i + j - 1)
        i = j
    T2 = np.empty(N, dtype=np.float)
    # Note(kazeevn) +1 is due to Python using 0-based indexing
    # instead of 1-based in the AUC formula in the paper
    T2[J] = T + 1
    return T2


def _fastDeLong(predictions_sorted_transposed, label_1_count):
    """
    The fast version of DeLong's method for computing the covariance of
    unadjusted AUC.
    Args:
       predictions_sorted_transposed: a 2D numpy.array[n_classifiers, n_examples]
          sorted such as the examples with label "1" are first
    Returns:
       (AUC value, DeLong covariance)
    Reference:
     @article{sun2014fast,
       title={Fast Implementation of DeLong's Algorithm for
              Comparing the Areas Under Correlated Receiver Oerating Characteristic Curves},
       author={Xu Sun and Weichao Xu},
       journal={IEEE Signal Processing Letters},
       volume={21},
       number={11},
       pages={1389--1393},
       year={2014},
       publisher={IEEE}
     }
    """
    # Short variables are named as they are in the paper
    m = label_1_count
    n = predictions_sorted_transposed.shape[1] - m
    positive_examples = predictions_sorted_transposed[:, :m]
    negative_examples = predictions_sorted_transposed[:, m:]
    k = predictions_sorted_transposed.shape[0]

    tx = np.empty([k, m], dtype=np.float)
    ty = np.empty([k, n], dtype=np.float)
    tz = np.empty([k, m + n], dtype=np.float)
    for r in range(k):
        tx[r, :] = _compute_midrank(positive_examples[r, :])
        ty[r, :] = _compute_midrank(negative_examples[r, :])
        tz[r, :] = _compute_midrank(predictions_sorted_transposed[r, :])
    aucs = tz[:, :m].sum(axis=1) / m / n - float(m + 1.0) / 2.0 / n
    v01 = (tz[:, :m] - tx[:, :]) / n
    v10 = 1.0 - (tz[:, m:] - ty[:, :]) / m
    sx = np.cov(v01)
    sy = np.cov(v10)
    delongcov = sx / m + sy / n
    return aucs, delongcov


def _calc_pvalue(aucs, sigma):
    """Computes log(10) of p-values.
    Args:
       aucs: 1D array of AUCs
       sigma: AUC DeLong covariances
    Returns:
       log10(pvalue)
    """
    l = np.array([[1, -1]])
    z = np.abs(np.diff(aucs)) / np.sqrt(np.dot(np.dot(l, sigma), l.T))
    return np.log10(2) + stats.norm.logsf(z, loc=0, scale=1) / np.log(10)


def _compute_ground_truth_statistics(ground_truth, sample_weight=None):
    assert np.array_equal(np.unique(ground_truth), [0, 1])
    order = (-ground_truth).argsort()
    label_1_count = int(ground_truth.sum())
    if sample_weight is None:
        ordered_sample_weight = None
    else:
        ordered_sample_weight = sample_weight[order]

    return order, label_1_count, ordered_sample_weight


def delong_roc_variance(ground_truth, predictions):
    """
    Computes ROC AUC variance for a single set of predictions
    Args:
       ground_truth: np.array of 0 and 1
       predictions: np.array of floats of the probability of being class 1
    """
    sample_weight = None
    order, label_1_count, ordered_sample_weight = _compute_ground_truth_statistics(
        ground_truth, sample_weight)
    predictions_sorted_transposed = predictions[np.newaxis, order]
    aucs, delongcov = _fastDeLong(predictions_sorted_transposed, label_1_count, ordered_sample_weight)
    assert len(aucs) == 1, "There is a bug in the code, please forward this to the developers"
    return aucs[0], delongcov


def delong_roc_test(ground_truth, predictions_one, predictions_two):
    """
    Computes log(p-value) for hypothesis that two ROC AUCs are different
    Args:
       ground_truth: np.array of 0 and 1
       predictions_one: predictions of the first model,
          np.array of floats of the probability of being class 1
       predictions_two: predictions of the second model,
          np.array of floats of the probability of being class 1
    """
    sample_weight = None
    order, label_1_count, _ = _compute_ground_truth_statistics(ground_truth)
    predictions_sorted_transposed = np.vstack((predictions_one, predictions_two))[:, order]
    aucs, delongcov = _fastDeLong(predictions_sorted_transposed, label_1_count)
    return _calc_pvalue(aucs, delongcov)


def volcano_plot(df, coef_col, pval_col, count_col,
                 s=60, alpha=0.5,
                 cutoff=None, cutoff_col=None, cutoff_label=None,
                 annot=None, annot_up_down=False, annot_size=10,
                 xlabel='log$_2$ fold change', ylabel='$-$log$_{10}$ pval',
                 ax=None):

    pval_col[pval_col] = -np.log10(df[pval_col])
    ax = sns.scatterplot(data=df, x=coef_col, y=pval_col, hue=count_col, s=s, alpha=alpha, ax=ax)

    if annot is not None:
        idx = []
        for direction in [1, -1] if annot_up_down else [None]:
            if direction is not None:
                annot_df = df.loc[df[coef_col] * direction > 0]
            else:
                annot_df = df
            idx.extend(annot_df.sort_values(pval_col, ascending=False).iloc[:annot].index.tolist())
        ax.annotate(df.loc[idx].index, (df.loc[idx, coef_col], df.loc[idx, pval_col]), fontsize=annot_size)

    if cutoff_col is not None and cutoff is not None and (df[cutoff_col] < cutoff).any():
        cutoff_on_y_axis = df.loc[(df[cutoff_col] - cutoff).abs.idxmin(), pval_col]
        ax.axhline(cutoff_on_y_axis, color='gray')
        ax.text(ax.get_xlim()[1] * 1.1, cutoff_on_y_axis,
                '{} {}'.format(cutoff_label, cutoff), fontsize=annot_size, va='center', ha='left', backgroundcolor='w')

    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    sns.despine()
    return ax


def replace_chars_for_LM(s, replace_chars=OLS_REPLACE_CHARS, replace_with='_'):
    for c in replace_chars:
        if isinstance(s, str):
            s = s.replace(c, replace_with)
        elif isinstance(s, (pd.Series, pd.Index)):
            s = s.str.replace(c, replace_with)
        elif is_iterable(s):
            s = pd.Series(s).str.replace(c, replace_with).tolist()
    return s


def drop_na_for_LM(X_df, Y_df, verbose=True):
    # remove samples where none of the Y variables was measured
    # remove smaples where any of the X was not measured
    assert X_df.index.equals(Y_df.index)
    if verbose:
        print('Removing {} any_null samples in X_df'.format(X_df.isnull().any(axis=1).sum()))
        print('Removing {} all_null samples in Y_df'.format(Y_df.isnull().all(axis=1).sum()))
    Y_df = Y_df.loc[~Y_df.isnull().all(axis=1) & ~X_df.isnull().any(axis=1)]
    X_df = X_df.loc[Y_df.index]
    return X_df, Y_df


def _filter_LM_coefs(coef, do_not_correct=None, filter_intercept=False):

    def _fullmatch(c):
        return re.fullmatch(c, coef) \
               or re.fullmatch(r'{}\[.*\]'.format(c), coef) \
               or re.fullmatch(r'C\({}\)\[.*\]'.format(c), coef)

    if filter_intercept and coef == 'Intercept':
        return False
    elif do_not_correct and any([_fullmatch(c) for c in do_not_correct]):
        return False
    else:
        return True


def _encode_coef(c, is_categorical, mark_categorical=True, add_categorical_suffix=False):
    encoded = ('C({}){}' if mark_categorical and is_categorical else '{}').format(
        replace_chars_for_LM(c), '\[.*\]' if add_categorical_suffix else '')
    if add_categorical_suffix:
        encoded = encoded.replace('(', '\(').replace(')', '\)')
    return encoded


def _rename_LM_coefs(coef, style='python'):
    if style == 'R':
        return 'Intercept' if coef == '(Intercept)' else re.sub(
            r'factor\(([A-Za-z0-9\_]+)\)([A-Za-z0-9\_]+)', r'\g<1>[T.\g<2>]', coef)
    elif style == 'python':
        return re.sub(r'C\(([A-Za-z0-9\_]+)\)([A-Za-z0-9\_\[\]\.]+)', r'\g<1>\g<2>', coef)
    else:
        raise ValueError


def _rename_python_coefs(coef):
    return _rename_LM_coefs(coef, style='python')


def _rename_R_coefs(coef):
    return _rename_LM_coefs(coef, style='R')


def correct_with_LM(model, y, do_not_correct=None, do_not_correct_intercept=True, verbose=True):
    _data = model.model.data
    mask = np.asarray([_filter_LM_coefs(coef, do_not_correct, filter_intercept=do_not_correct_intercept) for coef in _data.xnames])
    coefs_to_correct = np.asarray(_data.xnames)[mask]
    if verbose:
        print('\nCorrecting:\n{}'.format(coefs_to_correct))
        print('\nNOT correcting:\n{}\n'.format(np.asarray(_data.xnames)[~mask]))
    corrected_y = y - np.dot(_data.exog[:, mask], model.params.loc[coefs_to_correct])
    return corrected_y


def _get_dummy_LM_col(coef, data):
    if coef in data:
        x = data[coef].values
    else:
        c1, c2 = coef.split('[')
        c2 = c2[2:-1]
        x = pd.get_dummies(data, drop_first=True)['{}_{}'.format(c1, c2)].values
    return x


def _sum_group_var(var_df, var_groups):
    for g in var_groups:
        idx = var_df.index[var_df.index.str.startswith(g)]
        if len(idx) != 0:
            g_sum = var_df.loc[idx].sum()
            var_df.drop(idx, inplace=True)
            var_df.loc[g] = g_sum
    return var_df


def get_var_expl_for_fixeff_of_lmm(coefs, data, var_groups=None):
    var_df = pd.DataFrame()
    _var_groups = defaultdict(lambda: [])
    _expl_sum = 0
    for coef in coefs.index:
        if coef != 'Intercept':
            _expl = _get_dummy_LM_col(coef, data) * coefs.loc[coef]
            _var = np.var(_expl, ddof=0)
            var_df = var_df.append(pd.Series([_var], index=['Var'], name=coef))
            _expl_sum += _expl
    var_df = var_df['Var'] / var_df['Var'].sum() * np.var(_expl_sum, ddof=0)
    if var_groups:
        var_df = _sum_group_var(var_df, var_groups)
    return pd.DataFrame(var_df.rename('Var'))


def get_var_expl_for_coefs(model, var_groups=None, anova_typ=1):
    anova_df = anova_lm(model, typ=anova_typ)
    var_df = anova_df['sum_sq'] / anova_df['sum_sq'].sum()
    if var_groups:
        var_df = _sum_group_var(var_df, var_groups)
    return pd.DataFrame(var_df.rename('FVE'))


def fit_linear_model(X_df, Y_df, design=None, lmm_groups=None, scale=False, reml=True, variance_partition=False,
                     contrasts=None, permutation=None,
                     do_not_correct=None, return_corrected_X=False, var_groups=None,
                     just_correction=False, random_state=None, verbose=1):
    np.random.seed(random_state)
    if isinstance(do_not_correct, str):
        do_not_correct = [do_not_correct]
    if contrasts:
        assert isinstance(contrasts, str) or is_iterable(contrasts)
        if isinstance(contrasts, str):
            contrasts = [[contrasts]]
        elif not is_iterable(contrasts[0]):
            contrasts = [contrasts]
    if verbose > 1:
        print('Running the following linear model:')
    assert X_df.index.equals(Y_df.index)

    if permutation:
        if verbose:
            print('Permuting Y_df with random_state={}'.format(permutation))
        _index = Y_df.index.values
        Y_df = Y_df.sample(frac=1, random_state=permutation)
        Y_df.index = _index

    _F_test_cols = ['p.value', 'resid.df', 'stat']
    X_df = X_df.copy()
    X_df.columns = replace_chars_for_LM(X_df.columns)
    corrected_df = pd.DataFrame(columns=Y_df.index) if return_corrected_X else None
    coefs_df, pvals_df, tvals_df, SEs_df, Ns_df, df_resid_df, var_df = pd.DataFrame(), pd.DataFrame(), pd.DataFrame(), pd.DataFrame(), pd.DataFrame(), pd.DataFrame(), pd.DataFrame()
    contrasts_df = []
    lm_models = {}
    _first = True
    for y in Y_df.columns:
        _y = replace_chars_for_LM(y)
        nulls = Y_df[y].isnull()
        _skip = False
        for col in X_df.columns:
            if X_df.loc[~nulls, col].drop_duplicates().shape[0] < 2:
                print('ERROR: {} is unique for target {}, skipping...'.format(col, y))
                _skip = True
        if _y in X_df.columns:
            _skip = True
        try:
            if not _skip:
                if verbose:
                    print('y:', y)
                # booleans to nullable integers
                _cast_bools_and_ints = {c: float for c in X_df.columns if is_bool_dtype(X_df[c].dtype) or is_integer_dtype(X_df[c].dtype)}
                if X_df.shape[1] != 0:
                    _X_df = X_df.loc[~nulls].astype(_cast_bools_and_ints)
                    assert not _X_df.isnull().any().any()
                    if scale:
                        assert all([is_string_dtype(_X_df[c].dtype) or is_float_dtype(_X_df[c].dtype) for c in _X_df.columns]), _X_df.dtypes
                        _scalable_cols = [c for c in _X_df.columns if is_float_dtype(_X_df[c].dtype) and set(_X_df[c]) != set([0, 1])]
                        _X_df.loc[:, _scalable_cols] = StandardScaler().fit_transform(_X_df.loc[:, _scalable_cols].values)
                    data = pd.concat([_X_df, Y_df.loc[~nulls, y].rename(_y)], axis=1)
                else:
                    data = pd.DataFrame(Y_df.loc[~nulls, y].rename(_y))
                if lmm_groups is None:
                    formula = '{} ~ {}'.format(_y, design)
                    if (verbose == 1 and _first) or verbose > 1:
                        print('~ {}'.format(formula.split('~')[1].strip()))
                    lm_models[y] = ols(formula=formula, data=data).fit()

                    if contrasts:
                        _xnames = np.asarray([_rename_python_coefs(c) for c in lm_models[y].model.data.xnames])
                        _contrasts_df = pd.DataFrame(columns=_F_test_cols)
                        for contrast in contrasts:
                            encoded_contrast = np.asarray([np.asarray([re.match('{}(\[.*\])?'.format(c), coef) is not None for coef in _xnames]) for c in contrast])
                            contrast_name = '__'.join([coef for c in encoded_contrast for coef in _xnames[c.astype(bool)]])
                            F_test = lm_models[y].f_test(r_matrix=encoded_contrast.astype(float))
                            _contrasts_df = _contrasts_df.append(
                                pd.Series([F_test.pvalue, F_test.df_denom, F_test.fvalue[0][0]],
                                          index=_F_test_cols, name=contrast_name))

                    _coefs = lm_models[y].params.rename(_rename_python_coefs)
                    _pvals = lm_models[y].pvalues.rename(_rename_python_coefs)
                    _tvals = lm_models[y].tvalues.rename(_rename_python_coefs)
                    _SEs = lm_models[y].bse.rename(_rename_python_coefs)
                    _Ns = pd.Series([len(data)] * len(_coefs), index=_coefs.index)
                    _df_resid = pd.Series(
                        [getattr(lm_models[y], 'df_resid_inference', lm_models[y].df_resid)] * len(_coefs),
                        index=_coefs.index)

                    if variance_partition:
                        _var_df = get_var_expl_for_coefs(model=lm_models[y], var_groups=var_groups, anova_typ=1)
                        _var_df['target'] = y
                        _var_df.index.name = 'contrast'
                        _var_df = _var_df.reset_index().set_index(['target', 'contrast'])

                    if return_corrected_X:
                        corrected_y = correct_with_LM(model=lm_models[y], y=Y_df.loc[~nulls, y],
                                                      do_not_correct=do_not_correct,
                                                      do_not_correct_intercept=False,
                                                      verbose=(verbose == 1 and _first) or verbose > 1)
                else:
                    _groups = replace_chars_for_LM(lmm_groups)
                    formula = '{} ~ {} + {}'.format(_y, design.replace('C(', 'factor('), ' + '.join(['(1 | {})'.format(g) for g in _groups]))
                    if (verbose == 1 and _first) or verbose > 1:
                        print('~ {}'.format(formula.split('~')[1].strip()))
                    if not just_correction:
                        import pymer4
                        lm_models[y] = pymer4.Lmer(formula, data=data)
                        lm_models[y].fit(REML=reml, summarize=False)
                        assert lm_models[y].design_matrix.columns.equals(lm_models[y].coefs.index)
                        _coefs_df = lm_models[y].coefs.rename(_rename_R_coefs)

                        if contrasts:
                            encoded_contrasts, contrast_names = [], []
                            for contrast in contrasts:
                                encoded_contrasts.append(np.asarray([np.asarray([re.match('{}(\[.*\])?'.format(c), coef) is not None for coef in _coefs_df.index]) for c in contrast]).astype(float))
                                contrast_names.append('__'.join([coef for c in encoded_contrasts[-1] for coef in _coefs_df.index[c.astype(bool)]]))
                            _contrasts_df = lm_models[y].contest(contrast=encoded_contrasts)
                            _contrasts_df.index = contrast_names
                            _contrasts_df = _contrasts_df.rename({'DenDF': 'resid.df', 'F value': 'stat', 'Pr(>F)': 'p.value'}, axis=1)
                            _contrasts_df = _contrasts_df[_F_test_cols]

                        # _check_model = ols(formula='{} ~ {}'.format(_y, design), data=data).fit()
                        # assert set(_coefs_df.index) == set(_check_model.params.rename(_rename_python_coefs).index)

                        _coefs = _coefs_df['Estimate']
                        _pvals = _coefs_df['P-val']
                        _tvals = _coefs_df['T-stat']
                        _SEs = _coefs_df['SE']
                        _Ns = pd.Series([len(data)] * len(_coefs), index=_coefs.index)
                        _df_resid = _coefs_df['DF']

                        if variance_partition:
                            _var_df = lm_models[y].ranef_var.copy().drop(['Name', 'Std'], axis=1)
                            if len(_coefs.drop('Intercept')) != 0:
                                _var_df = pd.concat([
                                    _var_df, get_var_expl_for_fixeff_of_lmm(coefs=_coefs, data=data, var_groups=var_groups)
                                ])
                            _var_df['FVE'] = _var_df['Var'] / _var_df['Var'].sum()
                            _var_df['target'] = y
                            _var_df.index.name = 'contrast'
                            _var_df = _var_df.reset_index().set_index(['target', 'contrast'])
                    else:
                        assert not contrasts
                        _coefs = pd.Series(dtype=float)
                        _pvals = pd.Series(dtype=float)
                        _tvals = pd.Series(dtype=float)
                        _SEs = pd.Series(dtype=float)
                        _Ns = pd.Series(dtype=float)
                        _df_resid = pd.Series(dtype=float)
                        _var_df = pd.DataFrame()

                    if return_corrected_X:
                        assert len(_groups) == 1
                        _py_model = mixedlm(formula='{} ~ {}'.format(_y, design), data=data, groups=data[_groups[0]])
                        try:
                            _py_model = _py_model.fit(reml=reml)
                        except ConvergenceWarning as w:
                            if just_correction:
                                print('SKIPPING:', w)
                                _py_model = None

                        if _py_model is not None:
                            if not just_correction:
                                assert set(_coefs_df.index) == set(_py_model.params.drop(['Group Var']).rename(_rename_python_coefs).index)
                                assert ((_coefs_df['Estimate'] - _py_model.params.rename(_rename_python_coefs).loc[_coefs_df.index]).abs() < 1e-3).all(), (_coefs_df['Estimate'] - _py_model.params.rename(_rename_python_coefs).loc[_coefs_df.index]).abs().sort_values()
                                assert (_coefs_df[['P-val']].corrwith(_py_model.pvalues.rename(_rename_python_coefs).loc[_coefs_df.index]) > 0.98).all()

                            corrected_y = correct_with_LM(model=_py_model, y=Y_df.loc[~nulls, y],
                                                          do_not_correct=do_not_correct,
                                                          do_not_correct_intercept=False,
                                                          verbose=(verbose == 1 and _first) or verbose > 1)
                            if not any([_groups[0] in dnc for dnc in do_not_correct]):
                                _random_intercepts = pd.Series(
                                    [_py_model.random_effects[g].loc['Group'] for g in sorted(_py_model.random_effects.keys())],
                                    index=sorted(_py_model.random_effects.keys())
                                )
                                if (verbose == 1 and _first) or verbose > 1:
                                    with np.printoptions(threshold=5):
                                        print('Also correcting random intercepts:\n{}\n'.format(_random_intercepts.index.values))
                                corrected_y -= pd.get_dummies(data[_groups[0]]).dot(_random_intercepts)
                        else:
                            corrected_y = None

                coefs_df = coefs_df.append(_coefs.rename(y))
                pvals_df = pvals_df.append(_pvals.rename(y))
                tvals_df = tvals_df.append(_tvals.rename(y))
                SEs_df = SEs_df.append(_SEs.rename(y))
                Ns_df = Ns_df.append(_Ns.rename(y))
                df_resid_df = df_resid_df.append(_df_resid.rename(y))
                if variance_partition:
                    var_df = pd.concat([var_df, _var_df])
                if contrasts:
                    _contrasts_df.index.name = 'contrast'
                    _contrasts_df['target'] = y
                    _contrasts_df = _contrasts_df.reset_index().set_index(['target', 'contrast'])
                    contrasts_df.append(_contrasts_df)
                if return_corrected_X and corrected_y is not None:
                    corrected_df = corrected_df.append(corrected_y)
            else:
                coefs_df.loc[y, :] = np.nan
                pvals_df.loc[y, :] = np.nan
                tvals_df.loc[y, :] = np.nan
                SEs_df.loc[y, :] = np.nan
                Ns_df.loc[y, :] = np.nan
                df_resid_df.loc[y, :] = np.nan
                if return_corrected_X:
                    corrected_df.loc[y, :] = np.nan
        except Exception:
            print('ERROR with', y)
            traceback.print_exception(*sys.exc_info())
            coefs_df.loc[y, :] = np.nan
            pvals_df.loc[y, :] = np.nan
            tvals_df.loc[y, :] = np.nan
            SEs_df.loc[y, :] = np.nan
            Ns_df.loc[y, :] = np.nan
            df_resid_df.loc[y, :] = np.nan
            if return_corrected_X:
                corrected_df.loc[y, :] = np.nan
            print('End of ERROR with', y)
        _first = False
    assert not return_corrected_X or (corrected_df.T.shape[0] == Y_df.shape[0] and corrected_df.T.shape[1] <= Y_df.shape[1])
    assert coefs_df.columns.equals(tvals_df.columns) and coefs_df.columns.equals(SEs_df.columns) and coefs_df.columns.equals(Ns_df.columns) and coefs_df.columns.equals(pvals_df.columns) and coefs_df.columns.equals(df_resid_df.columns)

    results_df = pd.DataFrame(index=pd.MultiIndex.from_arrays([[], []], names=['target', 'contrast']))
    if not just_correction:
        for df, name in [(coefs_df, 'Coef'), (pvals_df, 'p.value'), (tvals_df, 'stat'), (SEs_df, 'error'), (Ns_df, 'sample.size'), (df_resid_df, 'resid.df')]:
            s = df.stack()
            s.index.names = ['target', 'contrast']
            s.name = name
            results_df = pd.concat([results_df, s], axis=1)

    if contrasts:
        contrasts_df = pd.concat(contrasts_df, axis=0)

    _return = [lm_models, results_df.sort_index(axis=1)]
    if contrasts:
        _return.append(contrasts_df.sort_index(axis=1))
    if return_corrected_X:
        _return.append(corrected_df.T)
    _return.append(var_df)
    return _return


def deprecated_regress_out(X_df, Y_df, variables, do_not_correct_prefixes=[], verbose=True):
    variables = np.asarray(variables)
    lm = LinearRegression(fit_intercept=True)
    assert not X_df[variables].isnull().any().any()
    assert not Y_df.isnull().any().any()
    assert X_df.index.equals(Y_df.index)
    X_df = pd.get_dummies(X_df[variables], drop_first=True)
    lm.fit(X=X_df, y=Y_df)
    coefs_df = pd.DataFrame(data=lm.coef_, index=Y_df.columns, columns=X_df.columns)

    cov_mask = np.asarray([not any([col.startswith(skip_col) for skip_col in do_not_correct_prefixes]) \
                           for col in X_df.columns])
    if verbose:
        print('Linear model variables:', coefs_df.columns)
        print('Correcting for:', coefs_df.columns[cov_mask])
        print('NOT correcting for:', coefs_df.columns[~cov_mask])

    corrected_Y_df = Y_df - np.dot(X_df.loc[:, cov_mask], coefs_df.loc[:, cov_mask].T)

    return corrected_Y_df, X_df, coefs_df


def calc_norm_factors(counts_df, method):
    assert method in ['TMM', 'TMMwsp', 'RLE', 'upperquartile']
    from rpy2.robjects import numpy2ri, pandas2ri, r
    numpy2ri.activate()
    pandas2ri.activate()
    r.source('calcNormFactors.R')
    return r.calcNormFactors(counts_df, method=method)


# credit: Andre Rendeiro
# def r2pandas_df(r_df):
#     """Make :class:`pandas.DataFrame` from a ``R`` dataframe given by :class:`rpy`."""
#     df = pd.DataFrame(np.asarray(r_df)).T
#     df.columns = [str(x) for x in r_df.colnames]
#     df.index = [str(x) for x in r_df.rownames]
#     return df


# credit: Andre Rendeiro
def lola(bed_files, universe_file, output_files, databases, split_collections=False, cpus=8, verbose=False):
    from rpy2.robjects import numpy2ri, pandas2ri, r
    from rpy2.robjects.packages import importr
    numpy2ri.activate()
    pandas2ri.activate()
    importr('LOLA')

    if isinstance(bed_files, str):
        bed_files = [bed_files]
    if isinstance(output_files, str):
        output_files = [output_files]
    if isinstance(databases, str):
        databases = [databases]

    assert len(bed_files) == len(output_files)

    if verbose:
        print('Reading universe file {}'.format(universe_file))
    universe = r('LOLA::readBed')(universe_file)

    if verbose:
        print('Loading region set databases')
    _regionDB = r('LOLA::loadRegionDB')(np.array(databases))

    results = {}
    for output_file, bed_file in zip(output_files, bed_files):
        if verbose:
            print('Reading BED file {}'.format(bed_file))
        user_set = r('LOLA::readBed')(bed_file)

        if verbose:
            print('Running LOLA testing for file {}'.format(bed_file))
        lola_results = r('LOLA::runLOLA')(user_set, universe, _regionDB, cores=cpus)
        if not isinstance(lola_results, pd.DataFrame):
            lola_results = pd.DataFrame(np.asarray(lola_results))

        if verbose:
            print('Saving all results for file {}'.format(bed_file))
        lola_results.to_csv(output_file, index=False, sep='\t')

        if split_collections:
            for collection in lola_results['collection'].drop_duplicates():
                if verbose:
                    print('Saving results for collection {}'.format(collection))
                lola_results[lola_results['collection'] == collection].to_csv(
                    '{}.{}.tsv'.format(output_file[:-len('.tsv')] if output_file.endswith('.tsv') else output_file,
                                       collection),
                    index=False, sep='\t')

        results[output_file] = lola_results

    return results


def enrichment_plot(sorted_enr, top_n=10, score='P-value', term_col='Term', convert_to_neg_log10=True,
                    thresholds=None, thr_kwargs=dict(Down={}, Up={}),
                    kind='barh', stacked=True, color=dict(Down=BLUE, Up=RED), force_equal_xlim=False,
                    xlim=None, hack_xticks=False,
                    legend=True, legend_kwargs=dict(), fdr_behind=False, ax=None, **kwargs):
    # to make sure down is blue and red is up
    assert set(sorted_enr.keys()).issubset(['Down', 'Up']), set(sorted_enr.keys())
    if isinstance(color, dict):
        color = [color[d] for d in set(sorted_enr.keys())]

    _plot = lambda _df, _ax, _kwargs: _df.plot(ax=_ax, kind=kind, stacked=stacked, color=color, **_kwargs)
    _neg_log_str = '$-$log$_{10}$'
    _rename_score = {
        'pValueLog': '{} P-value'.format(_neg_log_str),
        'oddsRatio': 'Odds Ratio',
        'support': 'Overlap',
        'qValue': 'Adjusted P-value',
        'FDR q-val': 'Adjusted P-value',
        'NOM p-val': 'P-value',
        'NES': 'Normalized enrichment score',
        'maxRnk': 'Max Rank'
    }

    enr_df = pd.concat([sorted_enr[direction].reset_index(drop=True).head(top_n).rename(
        lambda x: '{} {}'.format(direction, x), axis=1) for direction in set(sorted_enr.keys())], axis=1)
    enr_df.index.name = 'rank'
    for direction in set(sorted_enr.keys()):
        _term_col = '{} {}'.format(direction, term_col)
        enr_df[_term_col] = enr_df[_term_col].fillna('')

    if convert_to_neg_log10:
        thresholds = thresholds.copy() if thresholds else None
        for direction in set(sorted_enr.keys()):
            _score_col = '{} {}'.format(direction, score)
            _mask = ~enr_df[_score_col].isnull() & (enr_df[_score_col] != 0)
            enr_df.loc[_mask, '{} {} {}'.format(direction, _neg_log_str, score)] = -np.log10(enr_df.loc[_mask, _score_col])
            if thresholds and thresholds[direction] != 0:
                thresholds[direction] = -np.log10(thresholds[direction])
        score = '{} {}'.format(_neg_log_str, score)
    if 'Down' in set(sorted_enr.keys()):
        enr_df.loc[:, 'Down {}'.format(score)] *= -1
    enr_df = enr_df.sort_index(ascending=False)

    _df = enr_df[['{} {}'.format(direction, score) for direction in set(sorted_enr.keys())]].rename(lambda x: x.split()[0], axis=1)

    if 'Down' in set(sorted_enr.keys()):
        _df.index = enr_df['Down {}'.format(term_col)]
        kwargs['legend'] = legend and fdr_behind
        lax = _plot(_df, ax, kwargs)
        lax.set_xlabel(_rename_score.get(score, score))
    else:
        lax = ax

    if 'Up' in set(sorted_enr.keys()):
        _df.index = enr_df['Up {}'.format(term_col)]
        kwargs['legend'] = legend and not fdr_behind
        rax = _plot(_df, lax.twinx() if lax else None, kwargs)
        rax.axvline(0, c='0.15', lw=mpl.rcParams['axes.linewidth'])
    else:
        rax = ax

    if force_equal_xlim:
        max_xlim = max(np.abs(rax.get_xlim()).max(), np.abs(lax.get_xlim()).max())
        max_xlim += 0.05 * max_xlim
        for _ax in [rax, lax]:
            _ax.set_xlim((-max_xlim, max_xlim))

    if xlim:
        for _ax in [rax, lax]:
            _ax.set_xlim(xlim)

    if hack_xticks:
        a, b = rax.get_xlim()
        for _ax in [rax, lax]:
            _ax.xaxis.set_ticks(np.arange(int(np.ceil(a * 10)), int(np.floor(b * 10)) + 1, 1) / 10)
            _ax.xaxis.set_major_formatter(mpl.ticker.FormatStrFormatter('%0.1f'))

    if thresholds:
        _ax = lax if fdr_behind else rax
        # on lax, the line will be behind the bars, on rax, the line will be in front of the bars
        if 'Down' in set(sorted_enr.keys()):
            _ax.axvline(-1 * thresholds['Down'], **thr_kwargs['Down'])
        if 'Up' in set(sorted_enr.keys()):
            _ax.axvline(thresholds['Up'], **thr_kwargs['Up'])

    _ax = lax if lax else rax
    _ax.set_xticklabels(['{:g}'.format(abs(tick)) for tick in _ax.get_xticks()])

    if legend:
        _ax = lax if lax and fdr_behind else rax
        if thresholds and ('label' in thr_kwargs['Down'] or 'label' in thr_kwargs['Up']):
            handles, labels = _ax.get_legend_handles_labels()
            _ax.legend(handles[1:] + handles[:1], labels[1:] + labels[:1], **legend_kwargs)
        else:
            _ax.legend(**legend_kwargs)

    return lax, rax


def setup_plotting(style='ticks', context='notebook', font_scale=1, font_size=None,
                   use_tex=False, rc=None):
    assert context in PLOT_CONTEXTS
    sns.set_color_codes()
    sns.set_style(style)

    if rc is None:
        rc = {}
    if 'font.sans-serif' not in rc:
        rc['font.sans-serif'] = ['Helvetica', 'Arial', 'Liberation Sans', 'sans-serif']
    if 'legend.title_fontsize' not in rc:
        rc['legend.title_fontsize'] = (11 if font_size is None else font_size) * PLOT_CONTEXTS[context] * font_scale
    if 'figure.titlesize' not in rc:
        rc['figure.titlesize'] = (12 if font_size is None else font_size) * PLOT_CONTEXTS[context] * font_scale
    rc['text.usetex'] = use_tex
    rc['svg.fonttype'] = 'none'

    if font_size is not None:
        for font_key in ['axes.labelsize', 'axes.titlesize', 'legend.fontsize',
                         'xtick.labelsize', 'ytick.labelsize', 'font.size']:
            rc[font_key] = font_size

    context_object = sns.plotting_context(context, font_scale=font_scale, rc=rc)
    mpl.rcParams.update(rc)
    mpl.rcParams.update(context_object)

    return context_object


def _treat_as_option_with_no_value(value: t.Any) -> bool:
    return value is True


def _option_should_be_skipped(value: t.Any) -> bool:
    return value is False or value is None


class ArgumentUnparser:

    """
    Modified from https://github.com/mbdevpl/argunparse
    For license see https://github.com/mbdevpl/argunparse/blob/master/LICENSE

    For performing reverse operation to what argparse.ArgumentParser does."""

    # pylint: disable=too-many-arguments
    def __init__(
            self, short_opt: str = '-', long_opt: str = '--', opt_value: str = ' ') -> None:

        assert isinstance(short_opt, str)
        assert isinstance(long_opt, str)
        assert isinstance(opt_value, str)

        self._short_opt = short_opt
        self._long_opt = long_opt
        self._opt_value = opt_value

    def unparse_arg(self, arg: t.Any) -> str:
        """Convert an object into a string that can be used as a command-line argument."""
        if isinstance(arg, Iterable) and not isinstance(arg, string_types):
            arg = ' '.join(map(str, arg))
        else:
            if not isinstance(arg, str):
                arg = str(arg)
            if ' ' in arg:
                arg = repr(arg)
        if not arg:
            arg = '""'
        return arg

    def unparse_args(self, arguments: t.Sequence[t.Any],
                     *, to_list: bool = False) -> t.Union[str, t.List[str]]:
        """Convert list to string of command-line args."""
        unparsed = []
        for arg in arguments:
            unparsed.append(self.unparse_arg(arg))
        if to_list:
            return unparsed
        return ' '.join(unparsed)

    def unparse_option(self, key: str, value: t.Any,
                       *, to_list: bool = False) -> t.Union[str, t.List[str]]:
        """Convert a key-value pair into a string that can be used as a command-line option."""
        if _option_should_be_skipped(value):
            return [] if to_list else ''
        unparsed_key = '{}{}'.format(self._long_opt if len(key) > 1 else self._short_opt, key)
        if not _treat_as_option_with_no_value(value):
            unparsed_value = self.unparse_arg(value)
        if to_list and (self._opt_value == ' ' or _treat_as_option_with_no_value(value)):
            if _treat_as_option_with_no_value(value):
                return [unparsed_key]
            return [unparsed_key, unparsed_value]
        if not _treat_as_option_with_no_value(value):
            unparsed_option = '{}{}{}'.format(unparsed_key, self._opt_value, unparsed_value)
        if to_list:
            return [unparsed_option]
        if _treat_as_option_with_no_value(value):
            return unparsed_key
        return unparsed_option

    def unparse_options(self, options: t.Mapping[str, t.Any],
                        *, to_list: bool = False) -> t.Union[str, t.List[str]]:
        """Convert dictionary to string of command-line args."""
        unparsed = []
        for key, value in options.items():
            if _option_should_be_skipped(value):
                continue
            unparsed_option = self.unparse_option(key, value, to_list=to_list)
            if to_list:
                unparsed += unparsed_option
            else:
                unparsed.append(unparsed_option)
        if to_list:
            return unparsed
        return ' '.join(unparsed)

    def unparse_options_and_args(self, options: t.Mapping[str, t.Any], arguments: t.Sequence[t.Any],
                                 *, to_list: bool = False) -> t.Union[str, t.List[str]]:
        """Convert dictionary and list to string of command-line args."""
        if options is None:
            unparsed_options = [] if to_list else ''
        else:
            unparsed_options = self.unparse_options(options, to_list=to_list)
        if arguments is None:
            unparsed_args = [] if to_list else ''
        else:
            unparsed_args = self.unparse_args(arguments, to_list=to_list)
        if to_list:
            return unparsed_options + unparsed_args
        unparsed = []
        if unparsed_options:
            unparsed.append(unparsed_options)
        if unparsed_args:
            unparsed.append(unparsed_args)
        return ' '.join(unparsed)

    def unparse_to_list(self, *args, **kwargs) -> list:
        """Unparse given args as command-line arguments and kwargs as command-line options.

        Output is a list of strings.

        This process is a reverse of what built-in argparse module does with parse_args() method.
        """
        return self.unparse_options_and_args(kwargs, args, to_list=True)

    def unparse(self, *args, **kwargs) -> str:
        """Unparse given args as command-line arguments and kwargs as command-line options.

        This process is a reverse of what built-in argparse module does with parse_args() method.
        """
        return self.unparse_options_and_args(kwargs, args)


def read_RData(fn):
    import rpy2.robjects as robjects
    from rpy2.robjects import pandas2ri
    return pandas2ri.rpy2py_dataframe(robjects.r['readRDS'](fn))


def compare_F_test(ssr_restr, ssr_full, df_restr, df_full):
    assert df_restr > df_full
    df_diff = df_restr - df_full
    assert df_diff > 0
    assert (ssr_full > 0).all()
    assert (ssr_restr > 0).all()
    assert df_full > 0
    f_value = (ssr_restr - ssr_full) / df_diff / ssr_full * df_full
    p_value = stats.f.sf(f_value, df_diff, df_full)
    return f_value, p_value


def compare_F_test_from_n(ssr_restr, ssr_full, n_samples, k_restr, k_full):
    assert k_restr < k_full
    return compare_F_test(ssr_restr=ssr_restr, ssr_full=ssr_full,
                          df_restr=(n_samples - k_restr), df_full=(n_samples - k_full))


def compare_F_test_with_multiple_tests(ssr_restr, ssr_full, df_restr, df_full, multiple_test='fdr_bh'):
    assert isinstance(ssr_restr, pd.Series) and isinstance(ssr_full, pd.Series)
    assert ssr_restr.index.equals(ssr_full.index)
    F, pval = compare_F_test(ssr_restr=ssr_restr, ssr_full=ssr_full, df_restr=df_restr, df_full=df_full)
    df = pd.concat([F.rename('F'), pd.Series(pval, index=F.index, name='pval')], axis=1)
    if multiple_test:
        padj = multipletests(df['pval'], method=multiple_test)[1]
        df = pd.concat([df, pd.Series(padj, index=df.index, name='padj')], axis=1)
    return df


class MaskedPCA(BaseEstimator, TransformerMixin):

    def __init__(self, mask=None, n_components=None, copy=True, whiten=False, svd_solver='auto', tol=0.0,
                 iterated_power='auto', random_state=None):
        # a boolean mask containing selected features for PCA
        self.mask = mask
        self.pca = PCA(n_components=n_components, copy=copy, whiten=whiten, svd_solver=svd_solver, tol=tol,
                       iterated_power=iterated_power, random_state=random_state)

    def fit(self, X):
        _mask = self.mask if self.mask is not None else slice(None)
        self.pca.fit(X[:, _mask])
        return self

    def transform(self, X):
        _mask = self.mask if self.mask is not None else slice(None)
        _pca_features = self.pca.transform(X[:, _mask])
        if self.mask is not None:
            return np.hstack([X[:, ~_mask], _pca_features])
        else:
            return _pca_features


def adj_R2(R2, n, p):
    return 1 - (1 - R2) * (n - 1) / (n - p - 1)


def hue_jointplot(data, x, y, hue=None, hue_order=None, colormap=None,
                    figsize=None, fig=None, legend=True, scatter_kws=None, legend_kws=None):
    # defaults
    if colormap is None:
        colormap = sns.color_palette()  # ['blue','orange']
    if figsize is None:
        figsize = (5, 5)
    if fig is None:
        fig = plt.figure(figsize=figsize)
    if scatter_kws is None:
        scatter_kws = dict(alpha=0.5, lw=1)
    if legend_kws is None:
        legend_kws = dict()

    # derived variables
    if hue is None:
        return "use normal sns.jointplot"
    hue_groups = data[hue].unique()
    assert set(hue_order) == set(hue_groups)
    hue_groups = np.asarray(hue_order)

    subdata = dict()
    colors = dict()

    active_colormap = colormap[0: len(hue_groups)]
    #     legend_mapping = []
    for hue_grp, color in zip(hue_groups, active_colormap):
        #         legend_entry = mpatches.Patch(color=color, label=hue_grp)
        #         legend_mapping.append(legend_entry)

        subdata[hue_grp] = data[data[hue] == hue_grp]
        colors[hue_grp] = color

    # canvas setup
    grid = matplotlib.gridspec.GridSpec(2, 2,
                             width_ratios=[4, 1],
                             height_ratios=[1, 4],
                             hspace=0, wspace=0
                             )
    ax_main = plt.subplot(grid[1, 0])
    ax_xhist = plt.subplot(grid[0, 0], sharex=ax_main)
    ax_yhist = plt.subplot(grid[1, 1])  # , sharey=ax_main)

    ## plotting

    # histplot x-axis
    for hue_grp in hue_groups:
        sns.distplot(subdata[hue_grp][x], color=colors[hue_grp]
                     , ax=ax_xhist)

    # histplot y-axis
    for hue_grp in hue_groups:
        sns.distplot(subdata[hue_grp][y], color=colors[hue_grp]
                     , ax=ax_yhist, vertical=True)

    # main scatterplot
    # note: must be after the histplots else ax_yhist messes up
    #     for hue_grp in hue_groups:
    #         sns.regplot(data = subdata[hue_grp], fit_reg=False,
    #                     x = x, y = y, ax = ax_main, color = colors[hue_grp]
    #                     , scatter_kws=scatter_kws
    #                    )

    sns.scatterplot(ax=ax_main, data=data, x=x, y=y, hue=data[hue].values,
                    hue_order=hue_groups, palette=[colors[hue_grp] for hue_grp in hue_groups], **scatter_kws)
    ax_main.legend(**legend_kws).set_visible(legend)

    # despine
    for myax in [ax_yhist, ax_xhist]:
        sns.despine(ax=myax, bottom=False, top=True, left=False, right=True
                    , trim=False)
        plt.setp(myax.get_xticklabels(), visible=False)
        plt.setp(myax.get_yticklabels(), visible=False)

    # topright
    #     ax_legend   = plt.subplot(grid[0,1])#, sharey=ax_main)
    #     plt.setp(ax_legend.get_xticklabels(), visible=False)
    #     plt.setp(ax_legend.get_yticklabels(), visible=False)

    sns.despine(ax=ax_xhist, left=True)
    ax_xhist.set_xlabel(None)
    sns.despine(ax=ax_yhist, bottom=True)
    ax_yhist.set_ylabel(None)
    #     sns.despine(ax=ax_legend, bottom=True, left=True)

    #     ax_legend.legend(handles=legend_mapping)
    return fig, grid


def clustermap_legend(ax, color_codes, **kwargs):
    for value in color_codes:
        ax.bar(0, 0, linewidth=0, color=color_codes[value], label=value)
    ax.legend(**kwargs)


def heatmap_pval_annot(pvals, rows=None, cols=None, thr=0.05, use_asterisk=False, fmt='.0e',
                       asterisks=[0.05, 0.01, 0.001, 0.0001]):
    if isinstance(pvals, pd.Series):
        assert rows is not None and cols is not None
        assert len(pvals.index.names) == 2
        if rows.name == pvals.index.names[0] and cols.name == pvals.index.names[1]:
            get_idx = lambda x1, x2: (x1, x2)
        elif rows.name == pvals.index.names[1] and cols.name == pvals.index.names[0]:
            get_idx = lambda x1, x2: (x2, x1)
        else:
            raise ValueError
        pvals = np.asarray(
            [pvals.loc[get_idx(row, col)] if (col, row) in pvals.index else np.nan for row in rows for col in cols]
        ).reshape((len(rows), len(cols)))
    else:
        assert rows is None and cols is None
    assert pvals.ndim == 2
    pvals = pd.DataFrame(pvals)
    if thr:
        pvals = pvals.where(pvals <= thr)

    if use_asterisk:
        _format = lambda x: ''.join(['*' if x <= astx else '' for astx in set(asterisks)])
    else:
        fmt = '{' + (fmt if ':' in fmt else ':{}'.format(fmt)) + '}'
        _format = lambda x: fmt.format(x)
    pvals = pvals.applymap(lambda x: _format(x) if x and ~np.isnan(x) else '')
    return pvals


def corr_heatmap_with_pvals(df_x, df_y=None, corr_method='pearson', do_not_test_diag=False, lower_tria=False,
                            padj_method='fdr_bh', annot='pval', pval_asterisks=False, annot_fmt='.0e',
                            cluster=False, apply_X_clusters_to_Y=False, apply_Y_clusters_to_X=False, optimal_ordering=False,
                            cell_width=0.65, cell_height=0.5, xlabel=None, ylabel=None, title=None, title_offset=None, annot_size=None,
                            fig_fn=None, annot_kws={}, cbar_kws={}, **kwargs):
    assert annot in ['pval', 'corr', None], annot
    assert corr_method in ['pearson', 'spearman']
    if 'label' not in cbar_kws:
        cbar_kws['label'] = '{}{} correlation'.format(corr_method[0].upper(), corr_method[1:].lower())
    if annot_size is not None:
        annot_kws['size'] = annot_size
    if 'cmap' not in kwargs:
        kwargs['cmap'] = PALETTE_CODES['D']
    if 'fmt' not in kwargs:
        kwargs['fmt'] = '' if annot == 'pval' else annot_fmt
    if 'center' not in kwargs:
        kwargs['center'] = 0
    if 'xticklabels' not in kwargs:
        kwargs['xticklabels'] = 1
    if 'yticklabels' not in kwargs:
        kwargs['yticklabels'] = 1

    if df_y is None:
        df_y = df_x
        _check = df_x.corr(method=corr_method)
    else:
        if isinstance(df_y, pd.Series):
            df_y = pd.DataFrame(df_y)
        _check = None

    if corr_method == 'pearson':
        from scipy.stats import pearsonr
        corr_pval = lambda x, y: pearsonr(x, y)[1]
    elif corr_method == 'spearman':
        from scipy.stats import spearmanr
        corr_pval = lambda x, y: spearmanr(x, y)[1]

    corr_df, pval_df = [], []
    for col in df_y.columns:
        corr_df.append(df_x.corrwith(df_y[col], method=corr_method).rename(col))
        pval_df.append(df_x.corrwith(df_y[col], method=corr_pval).rename(col))
    corr_df = pd.concat(corr_df, axis=1).T
    pval_df = pd.concat(pval_df, axis=1).T
    if do_not_test_diag:
        pval_df = pval_df.where(~np.eye(df_y.shape[1], df_x.shape[1]).astype(bool))
    if _check is not None:
        pd.testing.assert_frame_equal(corr_df, _check)
    
    df = pd.concat([corr_df.stack().rename('R'), pval_df.stack().rename('p.value')], axis=1)
    # df.index.names = [ylabel, xlabel]
    if padj_method:
        from misc import adjusted_pvals
        df['padj'] = adjusted_pvals(df['p.value'], method=padj_method)
        padj_df = df['padj'].unstack()
        df['annot'] = df['padj']
    else:
        padj_df = None
        df['annot'] = df['p.value']

    R_df = df['R'].unstack().loc[corr_df.index, corr_df.columns]
    P_df = df['annot'].unstack().loc[corr_df.index, corr_df.columns]
    R_df.index.name = ylabel
    R_df.columns.name = xlabel
    figsize = (corr_df.shape[1] * cell_width, corr_df.shape[0] * cell_height)
    if cluster:
        from scipy.cluster.hierarchy import linkage
        plot_func = sns.clustermap
        kwargs['figsize'] = figsize

        _linkage_kws = dict(method=kwargs['method'] if 'method' in kwargs else 'average',
                           metric=kwargs['metric'] if 'metric' in kwargs else 'euclidean',
                           optimal_ordering=optimal_ordering)
        if apply_Y_clusters_to_X or apply_X_clusters_to_Y:
            if apply_Y_clusters_to_X:
                row_Z = linkage(R_df, **_linkage_kws)
                col_Z = row_Z
            elif apply_X_clusters_to_Y:
                col_Z = linkage(R_df.T, **_linkage_kws)
                row_Z = col_Z
        elif optimal_ordering:
            row_Z = linkage(R_df, **_linkage_kws)
            col_Z = linkage(R_df.T, **_linkage_kws)
        else:
            row_Z, col_Z = None, None
        kwargs['row_linkage'] = row_Z
        kwargs['col_linkage'] = col_Z
    else:
        if 'square' not in kwargs:
            kwargs['square'] = True
        plt.subplots(1, 1, figsize=figsize)
        plot_func = sns.heatmap
    use_asterisk, asterisks = False, {}
    if pval_asterisks:
        use_asterisk = True
        if is_iterable(pval_asterisks):
            asterisks = dict(asterisks=pval_asterisks)

    plot = plot_func(
        R_df, annot=heatmap_pval_annot(P_df, use_asterisk=use_asterisk, fmt=annot_fmt, **asterisks) if annot == 'pval' else True if annot == 'corr' else annot,
        cbar_kws=cbar_kws, annot_kws=annot_kws, **kwargs
    )
    if cluster:
        plot.fig.suptitle(title, y=title_offset if title_offset else 1)
        if 'aspect' in cbar_kws or 'shrink' in cbar_kws:
            x0, y0, x1, y1 = plot.cax.get_position().get_points().flatten()
            width, height = x1 - x0, y1 - y0
            if 'shrink' in cbar_kws:
                width = width * cbar_kws['shrink']
                height = height * cbar_kws['shrink']
            if 'aspect' in cbar_kws:
                height = height * (cbar_kws['aspect'] / (height / width))
            plot.cax.set_position([x0, y0, width, height])
        if lower_tria:
            mask = np.triu(np.ones_like(R_df), k=1).astype(bool)
            values = plot.ax_heatmap.collections[0].get_array().reshape(R_df.shape)
            values[mask] = np.nan
            plot.ax_heatmap.collections[0].set_array(values.flatten())
    else:
        plot.set_title(title, pad=title_offset)
    if fig_fn:
        savefig(fig_fn)
    return plot, corr_df, pval_df, padj_df


def stripplot(*args, **kwargs):
    if 'random_state' in kwargs:
        np.random.seed(kwargs.pop('random_state'))
    return sns.stripplot(*args, **kwargs)


# def sctransform(adata):
#     """
#     Function to call scTransform from Python
#     """
#     import scanpy as sc
#     from scipy.sparse import issparse
#     import rpy2.robjects as ro
#     import anndata2ri
#
#     ro.r('library(Seurat)')
#     anndata2ri.activate()
#
#     if issparse(adata.X):
#         if not adata.X.has_sorted_indices:
#             adata.X.sort_indices()
#
#     for key in adata.layers:
#         if issparse(adata.layers[key]):
#             if not adata.layers[key].has_sorted_indices:
#                 adata.layers[key].sort_indices()
#
#     ro.globalenv['adata'] = adata
#
#     ro.r('seurat_obj = as.Seurat(adata, counts="X", data=NULL)')
#
#     ro.r('res = SCTransform(object=seurat_obj, return.only.var.genes=FALSE, do.correct.umi=FALSE)')
#
#     norm_X = ro.r('res@assays$SCT@scale.data').T
#
#     return norm_X


def fisher_enrichment(regions1, regions2, universe):
    '''
    Enrichment of regions1 in regions2
                     regions2  not_regions2
        regions1         a        b
        not_regions1     c        d
    '''
    assert regions1.isin(universe).all()
    assert regions2.isin(universe).all()
    assert regions1.is_unique and regions2.is_unique and universe.is_unique

    a = len(regions1.intersection(regions2))
    b = len(regions1.difference(regions2))
    c = len(regions2.difference(regions1))
    d = len(universe.difference(regions1).difference(regions2))
    assert a + b + c + d == len(universe)

    odds_ratio, pval = fisher_exact([[a, b], [c, d]], alternative='greater')
    return odds_ratio, pval, f'{a}/{len(regions2)}'


def _nan_if_not_int(x):
    try:
        int(x)
        return x
    except:
        return np.nan


def regions_to_bed(regions, fn=None, skip_missing=True):
    assert is_iterable(regions)
    if not isinstance(regions, pd.Series):
        regions = pd.Series(regions)
    regions = regions.str.split(':', expand=True)
    assert regions.shape[1] == 3
    assert regions.index.is_unique
    regions[3] = regions.index
    regions.columns = ['chr', 'start', 'end', 'region_id']

    if skip_missing:
        regions['start'] = regions['start'].map(_nan_if_not_int, na_action='ignore')
        regions['end'] = regions['end'].map(_nan_if_not_int, na_action='ignore')
        null_mask = regions['chr'].isnull() | regions['start'].isnull() | regions['end'].isnull()
        if null_mask.sum() != 0:
            print(f'WARNING: {null_mask.sum()} null regions in file {fn}')
        regions = regions.loc[~null_mask]

    if fn is not None:
        regions.to_csv(fn, sep='\t', header=None, index=None)

    return regions


def _liftover(orig_regions, tmp_dir, bin_exec, over_chain):
    fn_orig = os.path.join(tmp_dir, 'tmp_300BCG_orig')
    fn_lifted = os.path.join(tmp_dir, 'tmp_300BCG_lifted')
    fn_unmapped = os.path.join(tmp_dir, 'tmp_300BCG_unmapped')

    non_null_orig_regions = regions_to_bed(regions=orig_regions, fn=fn_orig)
    n_nulls = len(orig_regions) - len(non_null_orig_regions)
    try:
        subprocess.run([
            bin_exec, fn_orig, over_chain, fn_lifted, fn_unmapped,
        ], check=True)
    except subprocess.CalledProcessError as e:
        raise e

    lifted_regions = pd.read_csv(fn_lifted, sep='\t', header=None, index_col=3)
    lifted_regions.index.name = orig_regions.index.name
    assert lifted_regions.shape[1] == 3
    lifted_regions = lifted_regions[0].astype(str) + ':' + \
                     lifted_regions[1].astype(str) + ':' + \
                     lifted_regions[2].astype(str)

    try:
        unmapped_regions = pd.read_csv(fn_unmapped, sep='\t', header=None, comment='#', index_col=3)
        unmapped_regions.index.name = orig_regions.index.name
        assert unmapped_regions.shape[1] == 3
        unmapped_regions = unmapped_regions[0].astype(str) + ':' + \
                           unmapped_regions[1].astype(str) + ':' + \
                           unmapped_regions[2].astype(str)
    except:
        unmapped_regions = None

    return lifted_regions, unmapped_regions, n_nulls


def liftover(orig_regions, remove_tmp=True,
             bin_exec='tools/liftOver/liftOver', over_chain='tools/liftOver/hg19ToHg38.over.chain.gz'):

    if remove_tmp:
        with TemporaryDirectory() as tmp_dir:
            return _liftover(orig_regions=orig_regions, tmp_dir=tmp_dir, bin_exec=bin_exec, over_chain=over_chain)
    else:
        return _liftover(orig_regions=orig_regions, tmp_dir='.', bin_exec=bin_exec, over_chain=over_chain)


def fill_in_gene_names(df, ens_to_gene):
    for idx in df.index:
        ens_id = df.loc[idx, 'human_gene_ensembl']
        if ens_id in ens_to_gene.index:
            df.loc[idx, 'human_gene_name'] = ens_to_gene.loc[ens_id]
    return df


def fill_in_ensembl_ids(df, ens_to_gene, organism, case_sensitive=False):
    assert organism in ['hs', 'human', 'mm', 'mouse']
    if organism == 'hs':
        organism = 'human'
    if organism == 'mm':
        organism = 'mouse'
    for idx in df.loc[df[f'{organism}_gene_ensembl'].isnull()].index:
        gene = df.loc[idx, f'{organism}_gene_name']
        if case_sensitive:
            match_mask = ens_to_gene == gene
        else:
            match_mask = ens_to_gene.str.lower() == gene.lower()
        ens_ids = ens_to_gene.index[match_mask]
        if len(ens_ids) != 0:
            df.loc[idx, f'{organism}_gene_ensembl'] = ';'.join(ens_ids)
    return df


def fix_ensembl_ids(df, ens_to_gene, ens_to_synonyms, organism, raise_missing=False, raise_mismatch=False, fix=True, verbose=1, warning_prefix=''):

    assert organism in ['hs', 'human', 'mm', 'mouse']
    if organism == 'hs':
        organism = 'human'
    if organism == 'mm':
        organism = 'mouse'

    if fix:
        gene_to_ens = {g.lower(): e for e, g in ens_to_gene.items() if not pd.isnull(g)}

    error_count_unknown_id = 0
    error_count_unknown_name = 0
    fixed_count_unknown_id = 0
    fixed_count_unknown_name = 0
    for idx in df.loc[~df[f'{organism}_gene_ensembl'].isnull()].index:
        ens_id = df.loc[idx, f'{organism}_gene_ensembl']
        gene = df.loc[idx, f'{organism}_gene_name']

        if fix:
            if ens_id not in ens_to_gene.index and not pd.isnull(gene) and gene.lower() in gene_to_ens.keys():
                fixed_count_unknown_id += 1
                if verbose > 1:
                    print(f'{warning_prefix}FIXED {ens_id} --> {gene_to_ens[gene.lower()]}')
                ens_id = gene_to_ens[gene.lower()]
                df.loc[idx, f'{organism}_gene_ensembl'] = ens_id
            if ens_id in ens_to_gene.index and gene != ens_to_gene[ens_id]:
                if pd.isnull(ens_to_gene[ens_id]) and not pd.isnull(gene):
                    print(gene, '-->', ens_to_gene[ens_id])
                fixed_count_unknown_name += 1
                if verbose > 1:
                    print(f"{warning_prefix}FIXED {ens_id}: {gene} --> {ens_to_gene[ens_id]}")
                gene = ens_to_gene[ens_id]
                df.loc[idx, f'{organism}_gene_name'] = gene

        if ens_id not in ens_to_gene.index:
            error_count_unknown_id += 1
            if verbose > 1:
                print(f'{warning_prefix}{ens_id}: invalid/retired Ensembl ID{f" ({gene})" if not pd.isnull(gene) else ""}')
            if raise_missing:
                raise ValueError
        else:
            synonyms = ens_to_synonyms.loc[ens_id]
            if not is_iterable(synonyms):
                synonyms = pd.Series([synonyms])
            if synonyms.isnull().any():
                assert synonyms.isnull().all()

            error = False
            if pd.isnull(gene):
                if not synonyms.isnull().all():
                    error = True
            else:
                if synonyms.isnull().all() or gene.lower() not in synonyms.str.lower().values:
                    error = True

            if error:
                error_count_unknown_name += 1
                if verbose > 1:
                    print(f"{warning_prefix}{ens_id}: {gene} not in {synonyms.tolist()}")
                if raise_mismatch:
                    print(f"ERROR: {warning_prefix}{ens_id}: {gene} not in {synonyms.tolist()}")
                    raise ValueError

    if verbose:
        if fixed_count_unknown_id != 0:
            print(f"{warning_prefix}FIXED {fixed_count_unknown_id}/{len(df)} invalid/retired Ensembl IDs")
        if fixed_count_unknown_name != 0:
            print(f"{warning_prefix}FIXED {fixed_count_unknown_name}/{len(df)} incorrect gene names")
        if error_count_unknown_id != 0:
            print(f"{warning_prefix}{error_count_unknown_id}/{len(df)} invalid/retired Ensembl IDs")
        if error_count_unknown_name != 0:
            print(f"{warning_prefix}{error_count_unknown_name}/{len(df)} incorrect gene names")

    return df


def get_ens_to_gene_map(fn, merge_genes_with_synonyms=True, chromosomes='auto', with_synonyms=True):
    if chromosomes == 'auto':
        chromosomes = [str(i) for i in range(1, 23)] + ['X', 'Y', 'MT']
    ens_to_gene = pd.read_csv(fn, sep='\t', usecols=['Gene stable ID', 'Gene name', 'Chromosome/scaffold name'])
    ens_to_gene = ens_to_gene.drop_duplicates()
    ens_to_gene = ens_to_gene.loc[ens_to_gene['Chromosome/scaffold name'].isin(chromosomes)]
    ens_to_gene = ens_to_gene.set_index('Gene stable ID')['Gene name']

    if with_synonyms:
        ens_to_synonyms = pd.read_csv(fn, sep='\t', usecols=['Gene stable ID', 'Gene Synonym', 'Chromosome/scaffold name'])
        ens_to_synonyms = ens_to_synonyms.drop_duplicates()
        ens_to_synonyms = ens_to_synonyms.loc[ens_to_synonyms['Chromosome/scaffold name'].isin(chromosomes)]
        ens_to_synonyms = ens_to_synonyms.set_index('Gene stable ID')['Gene Synonym']
        ens_to_synonyms = ens_to_synonyms.loc[~ens_to_synonyms.isnull()]
        if merge_genes_with_synonyms:
            ens_to_synonyms = pd.concat([ens_to_gene, ens_to_synonyms])
    else:
        if merge_genes_with_synonyms:
            ens_to_synonyms = ens_to_gene.copy()
        else:
            ens_to_synonyms = None

    return ens_to_gene, ens_to_synonyms


def mouse_to_human(df, mm_to_hs_ens):
    df['human_gene_ensembl'] = np.nan
    df['human_gene_name'] = np.nan
    for idx in df.loc[~df['mouse_gene_ensembl'].isnull()].index:
        mm_id = df.loc[idx, 'mouse_gene_ensembl']
        if mm_id in mm_to_hs_ens.index:
            df.loc[idx, 'human_gene_ensembl'] = mm_to_hs_ens.loc[mm_id, 'Gene stable ID']
            df.loc[idx, 'human_gene_name'] = mm_to_hs_ens.loc[mm_id, 'Gene name']
    return df


def get_mouse_human_homologs_MGI(fn, check_overlap=False, mm_to_hs_ens=None):
    # http://www.informatics.jax.org/faq/ORTH_dload.shtml
    assert not check_overlap or mm_to_hs_ens is not None
    df = pd.read_csv(fn, sep='\t')
    human_df = df.loc[df['Common Organism Name'] == 'human'].set_index('DB Class Key')
    human_df = human_df[~human_df.index.duplicated()]
    mouse_df = df.loc[df['Common Organism Name'] == 'mouse, laboratory'].set_index('DB Class Key')
    mouse_df = mouse_df[~mouse_df.index.duplicated()]
    mouse_df.columns = ['mm: ' + c for c in mouse_df.columns]
    map_df = human_df.join(mouse_df)
    assert not map_df['Symbol'].isnull().any()
    map_df = map_df.loc[~map_df['mm: Symbol'].isnull()]
    map_df = map_df.loc[~map_df['Symbol'].duplicated(keep=False)]

    match = 0
    no_match = []
    missing = []
    for idx in map_df.index:
        hs = map_df.loc[idx, 'Symbol'].lower()
        ms = map_df.loc[idx, 'mm: Symbol'].lower()
        ens_ms = mm_to_hs_ens.loc[mm_to_hs_ens['Gene name'].str.lower() == hs, 'Mouse gene name'].str.lower()
        if len(ens_ms) == 0:
            missing.append(hs)
        elif ms in ens_ms.values:
            match += 1
        else:
            no_match.append(hs)

    print(len(map_df), match, len(no_match), len(missing), (~mm_to_hs_ens['Mouse gene name'].str.lower().isin(map_df['mm: Symbol'].str.lower())).sum())
    # cca. (MGI: 17504, shared with Ensembl: 16096, disagree: 182, not in Ensembl: 1226, not in MGI: 2892)
    return map_df if not check_overlap else (map_df, match, no_match, missing)


def get_ens_mm_to_hs_map(fn, check=None):
    ens_mm_to_hs = pd.read_csv(fn, sep='\t')
    ens_mm_to_hs = ens_mm_to_hs.loc[~ens_mm_to_hs['Gene stable ID'].isnull() & ~ens_mm_to_hs['Mouse gene stable ID'].isnull()]
    ens_mm_to_hs = ens_mm_to_hs.loc[ens_mm_to_hs['Chromosome/scaffold name'].isin([str(i) for i in range(1, 23)] + ['X', 'Y', 'MT'])]
    ens_mm_to_hs = ens_mm_to_hs.loc[ens_mm_to_hs['Mouse chromosome/scaffold name'].isin([str(i) for i in range(1, 20)] + ['X', 'Y', 'MT'])]
    assert ens_mm_to_hs.set_index(['Gene stable ID', 'Mouse gene stable ID']).index.is_unique
    ens_mm_to_hs = ens_mm_to_hs.set_index('Mouse gene stable ID')
    if check is not None:
        assert check.loc[ens_mm_to_hs.index].equals(ens_mm_to_hs['Mouse gene name'])
    ens_mm_to_hs = ens_mm_to_hs.loc[~ens_mm_to_hs.index.duplicated(keep=False), ['Gene stable ID', 'Gene name']]
    return ens_mm_to_hs


def uropa_cmd(bed_fn, config_fn, gtf_fn, prefix, n_threads=None):
    cmd = "uropa --prefix {prefix} --input {config} --bed {bed} --gtf {gtf} {threads} --log {log}".format(
        prefix=prefix, config=config_fn, bed=bed_fn, gtf=gtf_fn,
        threads=f"--threads {n_threads}" if n_threads else "", log=f"{prefix}.log"
    )
    return cmd


def read_uropa_annotations(fn):
    annot_df = pd.read_csv(fn, sep='\t')
    annot_df = annot_df.set_index('peak_id')
    annot_df.loc[annot_df['feature'] == 'transcript', 'feat_type'] = \
        'transcript:' + annot_df.loc[annot_df['feature'] == 'transcript', 'transcript_type']
    annot_df.loc[annot_df['feature'] == 'gene', 'feat_type'] = \
        'gene:' + annot_df.loc[annot_df['feature'] == 'gene', 'gene_type']
    annot_df['length'] = annot_df['peak_end'] - annot_df['peak_start']
    annot_df.loc[annot_df['name'].isna(), 'name'] = 'NONE'

    annot_df = annot_df[
        ['peak_chr', 'peak_start', 'peak_end', 'length', 'feat_anchor', 'distance', 'relative_location', 'feat_type', 'gene_id', 'gene_name', 'name']
    ].rename({
        'peak_chr': 'chr', 'peak_start': 'start', 'peak_end': 'end', 'relative_location': 'location', 'name': 'characterization'
    }, axis=1)

    return annot_df


def run_uropa(regions, config_fn, gtf_fn, prefix, n_threads=None):
    if not isinstance(regions, pd.Series):
        regions = pd.Series(regions)

    with TemporaryDirectory() as tmp_dir:
        prefix = os.path.join(tmp_dir, prefix)
        bed_fn = f"{prefix}.bed"
        regions_to_bed(regions=regions, fn=bed_fn)

        cmd = uropa_cmd(bed_fn=bed_fn, config_fn=config_fn, gtf_fn=gtf_fn, prefix=prefix, n_threads=n_threads)
        try:
            subprocess.run([arg.strip() for arg in cmd.split()], check=True)
        except subprocess.CalledProcessError as e:
            raise e

        annot_df = read_uropa_annotations(f"{prefix}_finalhits.txt")

    assert regions.equals(
        annot_df['chr'].astype(str) + ':' + annot_df['start'].astype(str) + ':' + annot_df['end'].astype(str)
    )

    return annot_df


def simplify_genomic_locations(df, tss_proximal_label='TSS_proximal', gene_body_label='gene_body', distal_label='distal', intergenic_label='intergenic'):
    s = pd.Series([np.nan] * len(df), index=df.index, name='genomic_location', dtype=object)
    s.loc[df['characterization'].isin(TSS_PROXIMAL_GROUP)] = tss_proximal_label
    s.loc[df['characterization'].isin(GENE_BODY_GROUP)] = gene_body_label
    s.loc[df['characterization'].isin(DISTAL_GROUP)] = distal_label
    s.loc[df['characterization'].isin(INTERGENIC_GROUP)] = intergenic_label
    return s


def get_region_to_gene_mapping(peaks_df, region_filter, columns='gene_name', name=None):
    if isinstance(columns, str):
        columns = [columns]
    if name is None:
        name = \
            PROMOTER_MAP_COL if region_filter == TSS_PROXIMAL else \
        DISTAL_MAP_COL if region_filter == GENE_AND_DISTAL_10kb else None

    mapping = peaks_df[columns].copy()
    mask = peaks_to_genes(peaks_df, **PEAKS_TO_GENES[region_filter])
    mapping.loc[~mask, :] = np.nan
    if len(columns) == 1:
        mapping = mapping[columns[0]].rename(name)
    else:
        mapping.columns = [f"{name}:{c}" for c in mapping.columns]

    return mapping


def annotate_regions(regions, config_fn, gtf_fn, prefix, n_threads=None):
    annot_df = run_uropa(regions, config_fn, gtf_fn, prefix, n_threads=n_threads)
    extras = [simplify_genomic_locations(annot_df)]
    for region_filter in [TSS_PROXIMAL, GENE_AND_DISTAL_10kb]:
        extras.append(get_region_to_gene_mapping(annot_df, region_filter, columns=['gene_id', 'gene_name'],
                                                 name=region_filter))
    annot_df = pd.concat([annot_df] + extras, axis=1)
    return annot_df


def regions_to_genes(df, region_col, organism, config_fn, gtf_fn, prefix, n_threads=None, return_annot=False):
    assert organism in ['hs', 'mm']
    organism = 'human' if organism == 'hs' else 'mouse'

    annot_df = annotate_regions(
        regions=df.loc[~df[region_col].isnull(), region_col],
        config_fn=config_fn, gtf_fn=gtf_fn, prefix=prefix, n_threads=n_threads)
    annot_df = annot_df.reindex(df.index)

    assert annot_df[f"{TSS_PROXIMAL}:gene_id"].equals(
        annot_df[f"{GENE_AND_DISTAL_10kb}:gene_id"].where(annot_df['characterization'].isin(TSS_PROXIMAL_GROUP))
    )

    df[f"{organism}_gene_name"] = annot_df[f"{GENE_AND_DISTAL_10kb}:gene_name"]
    df[f"{organism}_gene_ensembl"] = annot_df[f"{GENE_AND_DISTAL_10kb}:gene_id"]
    df[f"{organism}_gene_ensembl"] = df[f"{organism}_gene_ensembl"].str.rsplit('.', n=1, expand=True)[0]
    assert df[f"{organism}_gene_name"].isnull().equals(df[f"{organism}_gene_ensembl"].isnull()), (df[f"{organism}_gene_name"].isnull().sum(), df[f"{organism}_gene_ensembl"].isnull().sum())

    df.insert(df.columns.tolist().index('hg38_region') + 1, 'promoter', False)
    df.loc[annot_df['characterization'].isin(TSS_PROXIMAL_GROUP), 'promoter'] = True
    df['promoter'] = df['promoter'].astype('boolean')
    df.loc[df[f"{organism}_gene_ensembl"].isnull(), 'promoter'] = np.nan

    return df if not return_annot else (df, annot_df)


def sclm(data_fn, celltype_col, results_dir, reference=None, interaction_model=False, subsample_celltypes=False):
    import anndata
    import statsmodels.formula.api as smf

    assert np.sum([reference is not None, interaction_model]) == 1
    assert np.sum([subsample_celltypes, interaction_model]) == 1
    MIN_TOTAL_CELLS = 10000
    adata = anndata.read_h5ad(data_fn)

    results = {}
    fn_prefix = f'DE_scLM_pearson_residuals.{"interaction" if interaction_model else reference}.{celltype_col}{".subsampled" if subsample_celltypes else ""}'
    celltypes = adata.obs[celltype_col].unique()
    celltypes = celltypes[~pd.isnull(celltypes)]

    if subsample_celltypes:
        only_these = []
        for celltype in celltypes:
            if (adata.obs[celltype_col] == celltype).sum() >= MIN_TOTAL_CELLS:
                only_these.append(celltype)
        celltypes = only_these

        df = adata.obs.loc[adata.obs[celltype_col].isin(celltypes)].copy()
        df[celltype_col] = df[celltype_col].astype(str)
        df = df.groupby([celltype_col, 'ts', 'ids'])['cells'].count().reset_index()
        n_cells = df.groupby(['ts', 'ids']).min()['cells']

    for celltype in celltypes:
        celltype_adata = adata[adata.obs[celltype_col] == celltype].copy()
        print(f'\nDIFFERENTIAL EXPRESSION: {celltype_col}, {"interaction" if interaction_model else reference}, {celltype} with {len(celltype_adata)} cells\n')

        if subsample_celltypes:
            idx = []
            for ts in celltype_adata.obs['ts'].unique():
                for ids in celltype_adata.obs['ids'].unique():
                    mask = (celltype_adata.obs['ts'] == ts) & (celltype_adata.obs['ids'] == ids)
                    assert mask.sum() >= n_cells.loc[(ts, ids)]
                    idx.append(np.random.RandomState(0).choice(
                        celltype_adata.obs_names[mask], size=n_cells.loc[(ts, ids)], replace=False))
            idx = np.concatenate(idx)
            celltype_adata = celltype_adata[celltype_adata.obs_names.isin(idx)]
            print('Subsampled:', len(celltype_adata))

        sys.stdout.flush()

        df = celltype_adata.obs[['ids', 'ts', 'time', 'stim']].copy()
        df['gene'] = np.nan
        df['gene'] = df['gene'].astype(float)

        assert df['ts'].isin(['T0_RPMI', 'T0_LPS', 'T3m_RPMI', 'T3m_LPS']).all()
        if interaction_model:
            model_formula = 'gene ~ ids + time * stim'
            df['time'] = df['time'].replace({'T0': f'AAAT0'})
            df['stim'] = df['stim'].replace({'RPMI': f'AAARPMI'})
        else:
            model_formula = 'gene ~ ids + ts'
            df['ts'] = df['ts'].replace({reference: f'AAA{reference}'})

        results[celltype] = {'coef': [], 'pval': [], 'tval': []}

        for g in range(celltype_adata.shape[1]):
            if g > 0 and g % 1000 == 0:
                print(g, end=', ')
                sys.stdout.flush()
            gene_X = celltype_adata.X[:, g]
            df['gene'] = np.asarray(gene_X.todense()).squeeze() if isinstance(gene_X, scipy.sparse.csr_matrix) else gene_X

            model = smf.ols(formula=model_formula, data=df).fit()
            results[celltype]['coef'].append(
                model.params.loc[model.params.index.str.contains('^ts\[|^time\[|^stim\[')].rename(celltype_adata.var_names[g]))
            results[celltype]['tval'].append(
                model.tvalues.loc[model.tvalues.index.str.contains('^ts\[|^time\[|^stim\[')].rename(celltype_adata.var_names[g]))
            results[celltype]['pval'].append(
                model.pvalues.loc[model.pvalues.index.str.contains('^ts\[|^time\[|^stim\[')].rename(celltype_adata.var_names[g]))

        print('all genes tested.')
        sys.stdout.flush()

        for k in results[celltype]:
            results[celltype][k] = pd.concat(results[celltype][k], axis=1).T
            results[celltype][k].to_csv(
                os.path.join(
                    results_dir,
                    f'{fn_prefix}.{celltype.replace(os.path.sep, "-")}.{k}.csv'
                )
            )

        del celltype_adata
        del df

    joblib.dump(results,
                os.path.join(results_dir, f'{fn_prefix}.results.pckl'),
                compress=True)


def sclm2(data_fn, celltype_col, results_dir, reference=None, interaction_model=False, subsample_celltypes=False, only_epi_genes=False, only_PBMCs=False, only_TRIM_genes=False):
    import anndata
    import statsmodels.formula.api as smf

    assert np.sum([reference is not None, interaction_model]) == 1
    assert np.sum([subsample_celltypes, interaction_model]) <= 1

    MIN_TOTAL_CELLS = 10000

    adata = anndata.read_h5ad(data_fn)
    print('adata:', adata.shape)

    if only_epi_genes:
        import misc
        peaks_df = misc.get_peak_annot(misc.PEAK_ANNOT_ALL_FN)
        epi_genes = peaks_df.loc[peaks_to_genes(peaks_df, **misc.PEAKS_TO_GENES[misc.GENE_AND_DISTAL_10kb]), 'gene_name'].values
        del peaks_df
        adata = adata[:, adata.var_names.isin(epi_genes)]
        print('adata:', adata.shape)

    if only_PBMCs:
        adata = adata[adata.obs['PBMC'] == 'PBMC']
        print('adata:', adata.shape)

    if only_TRIM_genes:
        TRIM_genes = pd.read_csv('results/scRNAseq/RNA_TRIM_genes.csv', header=None)[0].values
        adata = adata[:, adata.var_names.isin(TRIM_genes)]
        print('adata:', adata.shape)

    results = {}
    fn_prefix = f'DE_scLM_pearson_residuals.{"interaction" if interaction_model else reference}.{celltype_col}{".subsampled2" if subsample_celltypes else ""}'
    celltypes = adata.obs[celltype_col].unique()
    celltypes = celltypes[~pd.isnull(celltypes)]

    if subsample_celltypes:
        only_these = []
        for celltype in celltypes:
            if (adata.obs[celltype_col] == celltype).sum() >= MIN_TOTAL_CELLS:
                only_these.append(celltype)
        celltypes = only_these

        df = adata.obs.loc[adata.obs[celltype_col].isin(celltypes)].copy()
        df[celltype_col] = df[celltype_col].astype(str)
        df = df.groupby([celltype_col, 'ts'])['cells'].count().reset_index()
        n_cells = df.groupby(['ts']).min()['cells']

    for celltype in celltypes:
        if (adata.obs[celltype_col] == celltype).all():
            print('All just single cell type, memory saving :-)')
            celltype_adata = adata
        else:
            celltype_adata = adata[adata.obs[celltype_col] == celltype].copy()
        print(f'\nDIFFERENTIAL EXPRESSION: {celltype_col}, {"interaction" if interaction_model else reference}, {celltype} with {len(celltype_adata)} cells\n')

        if subsample_celltypes:
            idx = []
            for ts in celltype_adata.obs['ts'].unique():
                mask = (celltype_adata.obs['ts'] == ts)
                assert mask.sum() >= n_cells.loc[ts]
                idx.append(np.random.RandomState(0).choice(
                    celltype_adata.obs_names[mask], size=n_cells.loc[ts], replace=False))
            idx = np.concatenate(idx)
            celltype_adata = celltype_adata[celltype_adata.obs_names.isin(idx)]
            print('Subsampled:', len(celltype_adata))

        sys.stdout.flush()

        df = celltype_adata.obs[['ids', 'ts', 'time', 'stim']].copy()
        df['gene'] = np.nan
        df['gene'] = df['gene'].astype(float)

        assert df['ts'].isin(['T0_RPMI', 'T0_LPS', 'T3m_RPMI', 'T3m_LPS']).all()
        if interaction_model:
            model_formula = 'gene ~ ids + time * stim'
            df['time'] = df['time'].replace({'T0': f'AAAT0'})
            df['stim'] = df['stim'].replace({'RPMI': f'AAARPMI'})
        else:
            model_formula = 'gene ~ ids + ts'
            df['ts'] = df['ts'].replace({reference: f'AAA{reference}'})

        results[celltype] = {'coef': [], 'pval': [], 'tval': []}

        for g in range(celltype_adata.shape[1]):
            if g > 0 and g % 100 == 0:
                print(g, end=', ')
                sys.stdout.flush()
            gene_X = celltype_adata.X[:, g]
            df['gene'] = np.asarray(gene_X.todense()).squeeze() if isinstance(gene_X, scipy.sparse.csr_matrix) else gene_X

            model = smf.ols(formula=model_formula, data=df).fit()
            results[celltype]['coef'].append(
                model.params.loc[model.params.index.str.contains('^ts\[|^time\[|^stim\[')].rename(celltype_adata.var_names[g]))
            results[celltype]['tval'].append(
                model.tvalues.loc[model.tvalues.index.str.contains('^ts\[|^time\[|^stim\[')].rename(celltype_adata.var_names[g]))
            results[celltype]['pval'].append(
                model.pvalues.loc[model.pvalues.index.str.contains('^ts\[|^time\[|^stim\[')].rename(celltype_adata.var_names[g]))

        print('all genes tested.')
        sys.stdout.flush()

        for k in results[celltype]:
            results[celltype][k] = pd.concat(results[celltype][k], axis=1).T
            results[celltype][k].to_csv(
                os.path.join(
                    results_dir,
                    f'{fn_prefix}.{celltype.replace(os.path.sep, "-")}.{k}.csv'
                )
            )

        del celltype_adata
        del df

    joblib.dump(results,
                os.path.join(results_dir, f'{fn_prefix}.results.pckl'),
                compress=True)


# def bubble_heatmap():
#     # Plot as a bubble heatmap
#     _df = diff_enrich.loc[diff_enrich['description'].isin(terms)].copy()
#     _df['log10_padj_bh'] = -np.log10(_df['padj_bh'])
#     _df['log2_z_score'] = np.log2(_df['z_score'])  # z_score is actually
#     oddsRatio(misannotated)
#     _df['log2_z_score_dir'] = _df.apply(
#         lambda x:
#         (-1 if x['direction'] == 'down' else 1) * x['log2_z_score'],
#         axis=1)
#     _df['combined_score_dir'] = _df.apply(
#         lambda x:
#         (-1 if x['direction'] == 'down' else 1) * x['combined_score'],
#         axis=1)
#
#     # Assign x and y coordianates according to sorted comparison names
#     models = ['PBMC', 'Tcell', 'Bcell', 'Bcell-Naive', 'Bcell-Nateff',
#               'Bcell-Mem']
#     comps_x = [m for m in models if m in _df['model'].unique()]
#     comps_y = terms
#     gmap_x = {c: comps_x.index(c) for c in comps_x}
#     gmap_y = {c: comps_y.index(c) for c in comps_y}
#
#     _df['x'] = _df['model'].map(gmap_x)
#     _df['y'] = _df['description'].map(gmap_y)
#
#     # hue range
#     # hue_var = 'log2oddsRatio'
#     hue_var = 'log2_z_score_dir'
#     hue_min = au.round_down(_df[hue_var].min())
#     hue_max = au.round_up(_df[hue_var].max())
#     hue_range = max([hue_min, hue_max])
#
#     # Make figure perfectly square
#     h = 4
#     w = h / (len(comps_y) + 1) * (len(comps_x))
#     cbar_h = 0.25
#     sns.set_style('ticks')
#     fig, ax = plt.subplots(nrows=2,
#                            figsize=(w, h + cbar_h),
#                            gridspec_kw={'height_ratios': [h, cbar_h]})
#     sns.scatterplot(x='x',
#                     y='y',
#                     size='log10_padj_bh',
#                     sizes=(30, 300),
#                     size_norm=(1, 3),
#                     hue_norm=(-hue_range, hue_range),
#                     hue=hue_var,
#                     data=_df,
#                     palette='RdBu_r',
#                     edgecolor='black',
#                     linewidth=0.5,
#                     legend='brief',
#                     ax=ax[0])
#
#     ax[0].yaxis.tick_right()
#     ax[0].set_ylim(top=len(comps_y) - 0.5, bottom=-1)
#     ax[0].set_yticks([comps_y.index(c) for c in comps_y])
#     ax[0].set_yticklabels(comps_y)
#
#     ax[0].set_xlim(right=len(comps_x) - 0.5, left=-0.5)
#     ax[0].xaxis.tick_top()  # x axis on top
#     ax[0].xaxis.set_label_position('top')
#     ax[0].set_xticks([comps_x.index(c) for c in comps_x])
#     ax[0].set_xticklabels(comps_x, rotation=45, ha='left', va='bottom')
#
#     ax[0].invert_yaxis()
#     ax[0].tick_params(axis=u'both', which=u'both', length=0)
#     ax[0].set_ylabel('')
#     ax[0].set_xlabel('')
#
#     ax[0].legend(bbox_to_anchor=(0, -0.3), loc='upper left')
#
#     # Add a colorbar for log2oddsRatio
#     norm = plt.Normalize(-hue_range, hue_range)
#     sm = plt.cm.ScalarMappable(cmap="RdBu_r", norm=norm)
#     sm.set_array([])
#
#     # Remove the legend and add a colorbar
#     fig.colorbar(sm,
#                  orientation='horizontal',
#                  aspect=5,
#                  shrink=0.2,
#                  ticks=[hue_min, 0, hue_max],
#                  label='log2_z_score',
#                  cax=ax[1])
#
#     # Add spines as visual guides
#     for n in [comps_x.index(c) for c in comps_x]:
#         ax[0].axvline(n,
#                       linestyle='-',
#                       linewidth=0.5,
#                       color='black',
#                       alpha=0.3,
#                       zorder=0)
#     for n in [comps_y.index(c) for c in comps_y]:
#         ax[0].axhline(n,
#                       linestyle='-',
#                       linewidth=0.5,
#                       color='black',
#                       alpha=0.3,
#                       zorder=0)
#
#     # Save
#     sns.despine(fig, top=True, bottom=True, left=True, right=True)
#     # output_file = f'{output_dir}/{output_prefix}.{comparison}.cross_models.{threshold}.{step}.top_enrich.bubbleplot.svg'
#     # fig.savefig(output_file, bbox_inches='tight')