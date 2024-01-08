import os
import sys
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import numpy as np
import yaml
from scipy.cluster.hierarchy import linkage, fcluster
import misc
import bcg_utils as utils
from bcg_colors import *

_ = utils.setup_plotting(style='ticks', context='notebook', font_scale=1, rc=RC_PAPER)

config_fn = sys.argv[1]
with open(config_fn, 'r') as f:
    config = yaml.load(f, Loader=yaml.loader.FullLoader)
    model = config['model']
    celltype = config['celltype']
    results_dir = config['results_dir']

N_CLUSTERS = 5
FDR = 0.05

FIG_SIZE = (1.25, 1.25)
SAVE_FIG = True
fig_dir = os.path.split(misc.de_fn(celltype, model, results_dir=results_dir))[0]
clust_palette = [GREEN, YELLOW, PINK, CYAN, BROWN, RED, BLUE, ORANGE, PURPLE]
assert len(clust_palette) >= N_CLUSTERS

clusters = {}
for visit in ['V2', 'V3']:
    print(visit)
    extras = f'_DONOR.IC_DATE_2PI_SIN.{visit}_{visit}.DONOR.IC_DATE_2PI_COS'
    coefs_df = pd.read_csv(misc.de_fn(celltype, model, extras=extras, results_dir=results_dir), index_col=0)
    X_df = pd.read_csv(misc.de_fn(celltype, model, data='design', results_dir=results_dir), index_col=0)
    cos_coefs = coefs_df[misc.coef_col(f'{visit}.DONOR.IC_DATE_2PI_COS')].rename(f'{COEF_COL} (cosine)')
    sin_coefs = coefs_df[misc.coef_col(f'DONOR.IC_DATE_2PI_SIN.{visit}')].rename(f'{COEF_COL} (sine)')
    fitted_df = 0
    for c, coefs in [(f'{visit}.DONOR.IC_DATE_2PI_COS', cos_coefs),
                     (f'DONOR.IC_DATE_2PI_SIN.{visit}', sin_coefs)]:
        X = pd.DataFrame(data=np.broadcast_to(X_df[c], (len(coefs), len(X_df))).T, index=X_df.index, columns=coefs.index)
        fitted_df += X * coefs
    fitted_df = fitted_df.loc[fitted_df.index.str.contains(f'_{visit}_PBMC$')]

    df = misc.get_sample_annot()
    vacc_dates = df.loc[fitted_df.index, 'DONOR:IC_DATE_REAL'].sort_values()
    _first_vacc = 0.4
    _sort = vacc_dates % 1
    _sort.loc[_sort >= _first_vacc] -= 1
    vacc_dates = vacc_dates.iloc[np.argsort(_sort)]

    coef = f'SEASON.{visit}'
    de_df = misc.read_de(celltype, model, contrasts=coef, F_test=True, annot_fn=misc.PEAK_ANNOT_ALL_FN, results_dir=results_dir)
    print(f'{coef}: number of regions with <= FDR {FDR}: {(de_df[misc.padj_col(coef)] <= FDR).sum()}')

    # select the significant seasonal regions with a mask
    selected_seasonal_mask = de_df[misc.padj_col(coef)] < FDR
    fitted_df = fitted_df[de_df.index].loc[vacc_dates.index, selected_seasonal_mask]
    print(f'Fitted seasonal effects (FDR {FDR}): {fitted_df.shape}')

    # clustering
    cols_Z_fn = misc.de_fn(celltype, model, data=f'cols_Z_{visit}', ext='npz', gzipped=False, results_dir=results_dir)
    data_for_cols_Z_fn = misc.de_fn(celltype, model, data=f'data_for_cols_Z_{visit}', results_dir=results_dir)
    print('Calculating the linkage')
    cols_Z = linkage(fitted_df.T, method='average', metric='correlation', optimal_ordering=visit == 'V3')
    np.savez(cols_Z_fn, cols_Z=cols_Z)
    fitted_df.T.to_csv(data_for_cols_Z_fn)

    clusters[visit] = fcluster(cols_Z, N_CLUSTERS, criterion='maxclust')
    de_df[CLUSTER_COL] = np.nan
    de_df.loc[selected_seasonal_mask, CLUSTER_COL] = clusters[visit]
    de_df[CLUSTER_COL] = de_df[CLUSTER_COL].astype(pd.Int64Dtype())
    clustering_fn = misc.de_fn(celltype, model, data=f'clustering_{visit}', results_dir=results_dir)
    de_df.loc[~de_df[CLUSTER_COL].isnull(), [CLUSTER_COL]].astype(int).to_csv(clustering_fn)

    # number of regions per cluster
    fig, ax = plt.subplots(figsize=FIG_SIZE)
    sns.countplot(clusters[visit], palette=clust_palette)
    ax.set_xlabel(CLUSTER_COL)
    sns.despine()
    utils.savefig(os.path.join(fig_dir, f'atac_seasons_clusters_{visit}.pdf'), dpi=DPI)
    plt.close()

    # seasonal patterns
    peaks_df = misc.get_peak_annot()
    previous_xtick_dates = None
    fig, ax = plt.subplots(figsize=FIG_SIZE)
    for clust, color in zip(np.unique(clusters[visit]), clust_palette):
        _sort = np.argsort(fitted_df.loc[:, clusters[visit] == clust].mean(axis=1))[::-1]
        peak_month = misc.convert_partial_year(vacc_dates.iloc[_sort][0]).date().strftime("%b")
        print(f'Peak month for cluster {clust} ({(clusters[visit] == clust).sum()} regions): {peak_month}')
        ax.plot(vacc_dates.sort_values(), fitted_df.loc[vacc_dates.sort_values().index, clusters[visit] == clust].mean(axis=1), label=clust, c=color)
        ax.set_title('Mean of fitted seasonal effects')
        ax.set_ylabel('Seasonal contribution\nto the log2 fold change')
        xtick_dates = [misc.convert_partial_year(xtick).date() for xtick in ax.get_xticks()]
        if previous_xtick_dates is None:
            previous_xtick_dates = xtick_dates
        else:
            assert [str(d) for d in xtick_dates] == [str(d) for d in previous_xtick_dates]
    ax.set_xticklabels([d.strftime("%b %Y") for d in xtick_dates], rotation=30)
    ax.set_xlabel('Vaccination date')
    plt.legend(bbox_to_anchor=(1, 1), title=CLUSTER_COL)
    sns.despine()
    utils.savefig(os.path.join(fig_dir, f'atac_seasons_curves_{visit}.pdf'), dpi=DPI)
    plt.close()
