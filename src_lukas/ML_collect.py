import os
import warnings
from functools import wraps
import joblib

import matplotlib as mpl
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
from statsmodels.stats.multitest import multipletests
import re

import tables
import yaml

from misc import *
from bcg_utils import *


sns.set_color_codes()
sns.set_context('talk')
sns.set_style('ticks')
mpl.rcParams['legend.title_fontsize'] = 16
mpl.rcParams['figure.titlesize'] = 18

config_list_fn = sys.argv[1]
# simplify this
with open(os.path.join(ML_RESULTS_ROOT, config_list_fn) if not config_list_fn.startswith(ML_RESULTS_ROOT + os.path.sep) else config_list_fn) as f:
    config_filenames = list(map(lambda l: l.strip(), f.readlines()))

results_df = pd.DataFrame(columns=[
    'model', 'X_name', 'X', 'X_visits', 'y', 'y_visits', 'estimator',
    'features', 'DE_filter', 'select', 'scale', 'pca', 'poly', 'downsampling', 'binarize', 'seed',
    'scoring', 'alt_scoring', 'merge_folds', 'params', 'non_zero_coefs', 'coefs', 'test_size', 'class_ratio',
    'mean_cv_score', 'std_cv_score', 'mean_test_score', 'std_test_score', 'alt_mean_test_score', 'alt_std_test_score',
    'mean_train_score', 'std_train_score', 'union_test_score', 'alt_union_test_score',
    'FPR', 'TPR', 'precision', 'recall'])

for task_id in range(1, len(config_filenames) + 1):
    with open(config_filenames[task_id - 1], 'r') as f:
        config = yaml.load(f, Loader=yaml.loader.FullLoader)

    _dir = ML_RESULTS_ROOT if 'hack' not in config_list_fn else '{}_{}'.format(ML_RESULTS_ROOT, '_'.join(config_list_fn.split('_')[1:]))
    predictions_dir = os.path.join(_dir, get_data_name_from_config(config))
    print(predictions_dir)
    sys.stdout.flush()
    targets = get_numeric_targets_for_config(config, ML_RESULTS_ROOT)
    for target in targets:
        assert target.startswith(config['Y'].lstrip('^').rstrip('$')), (target, config['Y'].lstrip('^').rstrip('$'))
        results = evaluate_target(target, config, predictions_dir)
        if results is not None:
            results_df = df_append(results_df, name=None, **results)
        else:
            if target == 'SYSMEX:EO':
                pass
            else:
                print('Missing', get_filename_from_config(config, predictions_dir).format(data='{}.*'.format(target), ext='*'))

# usually several jobs are collecting different results and then they may try writing to the HDF file concurrently
success = False
for _ in range(5):
    try:
        with warnings.catch_warnings():
            warnings.simplefilter(action='ignore', category=tables.NaturalNameWarning)
            warnings.simplefilter(action='ignore', category=pd.errors.PerformanceWarning)
            results_df.to_hdf(os.path.join(ML_RESULTS_ROOT, '{}.ML.results.hdf'.format(PROJECT)), key=config_list_fn, complevel=9)
        success = True
    except:
        time.sleep(120)

# if could not write into the shared file, save in its own file
if not success:
    with warnings.catch_warnings():
        warnings.simplefilter(action='ignore', category=tables.NaturalNameWarning)
        warnings.simplefilter(action='ignore', category=pd.errors.PerformanceWarning)
        results_df.to_hdf(os.path.join(ML_RESULTS_ROOT, '{}.ML.results.{}.hdf'.format(PROJECT, config_list_fn)), key='results', complevel=9)
