#!/usr/bin/env python
# coding: utf-8

import matplotlib as mpl
mpl.use('Agg')
import tables

from pandas.api.types import is_numeric_dtype
from misc import *
from bcg_utils import *
from bcg_utils import _encode_coef

_ = setup_plotting(style='ticks', context='notebook')

SAVE_DATA = True
model = sys.argv[1] if len(sys.argv) > 1 else 'final_V1_corrected_combat'
print(model)

# _ATAC_COVS = ['LAB:BATCH', 'RUN:TSS_ENRICHMENT']
_SEX_AGE_BMI_CONTRA_TIME_ALCO = ['DONOR:SEX', 'DONOR:AGE', 'DONOR:BMI', 'DONOR:oralContraceptivesIncludingMen', 'SAMPLE:VISIT_TIME_REAL', 'SAMPLE:alcoholInLast24h']
_VISIT_SEX_AGE_BMI_CONTRA_TIME_ALCO = ['SAMPLE:VISIT'] + _SEX_AGE_BMI_CONTRA_TIME_ALCO

MODELS = {
    # 'corrected': {
    #     'ATAC': BLOOD + _VISIT_SEX_AGE_TIME + _ATAC_COVS,
    #     'CM': WHOLE_BLOOD + _VISIT_SEX_AGE_TIME + ['DONOR:BMI', 'DONOR:oralContraceptivesIncludingMen'],
    #     'CYTO': BLOOD + _VISIT_SEX_AGE_TIME,
    #     'PBMC_PERC': _VISIT_SEX_AGE_TIME,
    #     'PBMC_PER_ML': _VISIT_SEX_AGE_TIME,
    #     'method': None,
    #     'n_components': None,
    #     'visits': None,
    #     'get_data': get_norm_counts
    # },
    'blood_TSS_enr_corrected_combat': {
        'ATAC': BLOOD + ['SAMPLE:VISIT', 'SAMPLE:DONOR'] + ['RUN:TSS_ENRICHMENT'],
        'CM': None,
        'CYTO': None,
        'PBMC_PERC': None,
        'PBMC_PER_ML': None,
        'WB_PER_ML': None,
        'do_not_correct': ['SAMPLE:VISIT', 'SAMPLE:DONOR'],
        'method': None,
        'n_components': None,
        'visits': None,
        'get_data': get_batch_corrected_counts
    },
    'TSS_enr_corrected_combat': {
        'ATAC': BLOOD + ['SAMPLE:VISIT', 'SAMPLE:DONOR'] + ['RUN:TSS_ENRICHMENT'],
        'CM': None,
        'CYTO': None,
        'PBMC_PERC': None,
        'PBMC_PER_ML': None,
        'WB_PER_ML': None,
        'do_not_correct': BLOOD + ['SAMPLE:VISIT', 'SAMPLE:DONOR'],
        'method': None,
        'n_components': None,
        'visits': None,
        'get_data': get_batch_corrected_counts
    },
    'corrected_LMM_combat': {
        'do_not_return_without_convergence': False,
        'ATAC': BLOOD + _VISIT_SEX_AGE_BMI_CONTRA_TIME_ALCO + ['RUN:TSS_ENRICHMENT'],
        'CM': None,  # WHOLE_BLOOD + _VISIT_SEX_AGE_BMI_CONTRA_TIME_ALCO,
        'CYTO': None,  # BLOOD + _VISIT_SEX_AGE_BMI_CONTRA_TIME_ALCO,
        'PBMC_PERC': None,  # _VISIT_SEX_AGE_BMI_CONTRA_TIME_ALCO,
        'PBMC_PER_ML': None,  # _VISIT_SEX_AGE_BMI_CONTRA_TIME_ALCO,
        'WB_PER_ML': None,  # _VISIT_SEX_AGE_BMI_CONTRA_TIME_ALCO,
        'do_not_correct': ['SAMPLE_VISIT', 'SAMPLE_DONOR', 'Intercept'],
        'method': None,
        'n_components': None,
        'visits': None,
        'get_data': get_batch_corrected_counts
    },
    'corrected_converged_LMM_combat': {
        'do_not_return_without_convergence': True,
        'ATAC': BLOOD + _VISIT_SEX_AGE_BMI_CONTRA_TIME_ALCO + ['RUN:TSS_ENRICHMENT'],
        'CM': None,  # WHOLE_BLOOD + _VISIT_SEX_AGE_BMI_CONTRA_TIME_ALCO,
        'CYTO': None,  # BLOOD + _VISIT_SEX_AGE_BMI_CONTRA_TIME_ALCO,
        'PBMC_PERC': None,  # _VISIT_SEX_AGE_BMI_CONTRA_TIME_ALCO,
        'PBMC_PER_ML': None,  # _VISIT_SEX_AGE_BMI_CONTRA_TIME_ALCO,
        'WB_PER_ML': None,  # _VISIT_SEX_AGE_BMI_CONTRA_TIME_ALCO,
        'do_not_correct': ['SAMPLE_VISIT', 'SAMPLE_DONOR', 'Intercept'],
        'method': None,
        'n_components': None,
        'visits': None,
        'get_data': get_batch_corrected_counts
    },
    'final_corrected_combat': {
        'which': '.corrected',
        'ATAC': BLOOD + _VISIT_SEX_AGE_BMI_CONTRA_TIME_ALCO + ['RUN:TSS_ENRICHMENT'],
        'CM': WHOLE_BLOOD + _VISIT_SEX_AGE_BMI_CONTRA_TIME_ALCO,
        'CYTO': BLOOD + _VISIT_SEX_AGE_BMI_CONTRA_TIME_ALCO,
        'PBMC_PERC': _VISIT_SEX_AGE_BMI_CONTRA_TIME_ALCO,
        'PBMC_PER_ML': _VISIT_SEX_AGE_BMI_CONTRA_TIME_ALCO,
        'WB_PER_ML': _VISIT_SEX_AGE_BMI_CONTRA_TIME_ALCO,
        'do_not_correct': ['SAMPLE:VISIT_V2', 'SAMPLE:VISIT_V3'],
        'method': None,
        'n_components': None,
        'visits': None,
        'get_data': get_batch_corrected_counts
    },
    'final_V1_corrected_combat_alcoRem': {
        'which': '.corrected_alcoRem',
        'ATAC': BLOOD + _SEX_AGE_BMI_CONTRA_TIME_ALCO + ['RUN:TSS_ENRICHMENT'],
        'CM': WHOLE_BLOOD + _SEX_AGE_BMI_CONTRA_TIME_ALCO,
        'CYTO': BLOOD + _SEX_AGE_BMI_CONTRA_TIME_ALCO,
        'PBMC_PERC': _SEX_AGE_BMI_CONTRA_TIME_ALCO,
        'PBMC_PER_ML': _SEX_AGE_BMI_CONTRA_TIME_ALCO,
        'WB_PER_ML': _SEX_AGE_BMI_CONTRA_TIME_ALCO,
        'do_not_correct': [],
        'method': None,
        'n_components': None,
        'visits': ['V1'],
        'get_data': get_batch_corrected_counts
    },
    'final_V1_corrected_combat': {
        'which': '.corrected',
        'ATAC': BLOOD + _SEX_AGE_BMI_CONTRA_TIME_ALCO + ['RUN:TSS_ENRICHMENT'],
        'CM': WHOLE_BLOOD + _SEX_AGE_BMI_CONTRA_TIME_ALCO,
        'CYTO': BLOOD + _SEX_AGE_BMI_CONTRA_TIME_ALCO,
        'PBMC_PERC': _SEX_AGE_BMI_CONTRA_TIME_ALCO,
        'PBMC_PER_ML': _SEX_AGE_BMI_CONTRA_TIME_ALCO,
        'WB_PER_ML': _SEX_AGE_BMI_CONTRA_TIME_ALCO,
        'do_not_correct': [],
        'method': None,
        'n_components': None,
        'visits': ['V1'],
        'get_data': get_batch_corrected_counts
    },
    'final_corrected_combat_adjScores': {
        'which': '.corrected_adjScores',
        'ATAC': BLOOD + _VISIT_SEX_AGE_BMI_CONTRA_TIME_ALCO + ['RUN:TSS_ENRICHMENT'],
        'CM': WHOLE_BLOOD + _VISIT_SEX_AGE_BMI_CONTRA_TIME_ALCO,
        'CYTO': BLOOD + _VISIT_SEX_AGE_BMI_CONTRA_TIME_ALCO,
        'PBMC_PERC': _VISIT_SEX_AGE_BMI_CONTRA_TIME_ALCO,
        'PBMC_PER_ML': _VISIT_SEX_AGE_BMI_CONTRA_TIME_ALCO,
        'WB_PER_ML': _VISIT_SEX_AGE_BMI_CONTRA_TIME_ALCO,
        'do_not_correct': ['SAMPLE:VISIT_V2', 'SAMPLE:VISIT_V3'],
        'method': None,
        'n_components': None,
        'visits': None,
        'get_data': get_batch_corrected_counts
    },
    'corrected_combat_noS_adjScores': {
        'ATAC': BLOOD + _VISIT_SEX_AGE_BMI_CONTRA_TIME_ALCO + ['RUN:TSS_ENRICHMENT'],
        'CM': WHOLE_BLOOD + _VISIT_SEX_AGE_BMI_CONTRA_TIME_ALCO,
        'CYTO': BLOOD + _VISIT_SEX_AGE_BMI_CONTRA_TIME_ALCO,
        'PBMC_PERC': _VISIT_SEX_AGE_BMI_CONTRA_TIME_ALCO,
        'PBMC_PER_ML': _VISIT_SEX_AGE_BMI_CONTRA_TIME_ALCO,
        'WB_PER_ML': _VISIT_SEX_AGE_BMI_CONTRA_TIME_ALCO,
        'do_not_correct': ['SAMPLE:VISIT_V2', 'SAMPLE:VISIT_V3'],
        'method': None,
        'n_components': None,
        'visits': None,
        'get_data': get_batch_corrected_counts
    },
    'corrected_combat_noS': {
        'ATAC': BLOOD + _VISIT_SEX_AGE_BMI_CONTRA_TIME_ALCO + ['RUN:TSS_ENRICHMENT'] + ['SAMPLE:VISIT_DATE_2PI_COS', 'SAMPLE:VISIT_DATE_2PI_SIN'],
        'CM': WHOLE_BLOOD + _VISIT_SEX_AGE_BMI_CONTRA_TIME_ALCO + ['SAMPLE:VISIT_DATE_2PI_COS', 'SAMPLE:VISIT_DATE_2PI_SIN'],
        'CYTO': BLOOD + _VISIT_SEX_AGE_BMI_CONTRA_TIME_ALCO + ['SAMPLE:VISIT_DATE_2PI_COS', 'SAMPLE:VISIT_DATE_2PI_SIN'],
        'PBMC_PERC': _VISIT_SEX_AGE_BMI_CONTRA_TIME_ALCO + ['SAMPLE:VISIT_DATE_2PI_COS', 'SAMPLE:VISIT_DATE_2PI_SIN'],
        'PBMC_PER_ML': _VISIT_SEX_AGE_BMI_CONTRA_TIME_ALCO + ['SAMPLE:VISIT_DATE_2PI_COS', 'SAMPLE:VISIT_DATE_2PI_SIN'],
        'WB_PER_ML': _VISIT_SEX_AGE_BMI_CONTRA_TIME_ALCO + ['SAMPLE:VISIT_DATE_2PI_COS', 'SAMPLE:VISIT_DATE_2PI_SIN'],
        'do_not_correct': ['SAMPLE:VISIT_V2', 'SAMPLE:VISIT_V3'],
        'method': None,
        'n_components': None,
        'visits': None,
        'get_data': get_batch_corrected_counts
    },
    'corrected_combat_SX': {
        'ATAC': BLOOD + _VISIT_SEX_AGE_BMI_CONTRA_TIME_ALCO + ['RUN:TSS_ENRICHMENT'],
        'CM': WHOLE_BLOOD + _VISIT_SEX_AGE_BMI_CONTRA_TIME_ALCO,
        'CYTO': BLOOD + _VISIT_SEX_AGE_BMI_CONTRA_TIME_ALCO,
        'PBMC_PERC': _VISIT_SEX_AGE_BMI_CONTRA_TIME_ALCO,
        'PBMC_PER_ML': _VISIT_SEX_AGE_BMI_CONTRA_TIME_ALCO,
        'WB_PER_ML': _VISIT_SEX_AGE_BMI_CONTRA_TIME_ALCO,
        'do_not_correct': ['SAMPLE:VISIT_V2', 'SAMPLE:VISIT_V3'],
        'method': None,
        'n_components': None,
        'visits': None,
        'get_data': get_batch_corrected_counts
    },
    'PCA10_corrected_combat': {
        'ATAC': BLOOD + _VISIT_SEX_AGE_BMI_CONTRA_TIME_ALCO + ['RUN:TSS_ENRICHMENT'],
        'CM': WHOLE_BLOOD + _VISIT_SEX_AGE_BMI_CONTRA_TIME_ALCO,
        'CYTO': BLOOD + _VISIT_SEX_AGE_BMI_CONTRA_TIME_ALCO,
        'PBMC_PERC': _VISIT_SEX_AGE_BMI_CONTRA_TIME_ALCO,
        'PBMC_PER_ML': _VISIT_SEX_AGE_BMI_CONTRA_TIME_ALCO,
        'WB_PER_ML': _VISIT_SEX_AGE_BMI_CONTRA_TIME_ALCO,
        'do_not_correct': ['SAMPLE:VISIT_V2', 'SAMPLE:VISIT_V3'],
        'method': 'PCA',
        'n_components': 10,
        'visits': None,
        'get_data': get_batch_corrected_counts
    },
    'UMAP10_corrected_combat': {
        'ATAC': BLOOD + _VISIT_SEX_AGE_BMI_CONTRA_TIME_ALCO + ['RUN:TSS_ENRICHMENT'],
        'CM': WHOLE_BLOOD + _VISIT_SEX_AGE_BMI_CONTRA_TIME_ALCO,
        'CYTO': BLOOD + _VISIT_SEX_AGE_BMI_CONTRA_TIME_ALCO,
        'PBMC_PERC': _VISIT_SEX_AGE_BMI_CONTRA_TIME_ALCO,
        'PBMC_PER_ML': _VISIT_SEX_AGE_BMI_CONTRA_TIME_ALCO,
        'WB_PER_ML': _VISIT_SEX_AGE_BMI_CONTRA_TIME_ALCO,
        'do_not_correct': ['SAMPLE:VISIT_V2', 'SAMPLE:VISIT_V3'],
        'method': 'UMAP',
        'n_components': 10,
        'visits': None,
        'get_data': get_batch_corrected_counts
    },
    'PCA100_corrected_combat': {
        'ATAC': BLOOD + _VISIT_SEX_AGE_BMI_CONTRA_TIME_ALCO + ['RUN:TSS_ENRICHMENT'],
        'CM': WHOLE_BLOOD + _VISIT_SEX_AGE_BMI_CONTRA_TIME_ALCO,
        'CYTO': BLOOD + _VISIT_SEX_AGE_BMI_CONTRA_TIME_ALCO,
        'PBMC_PERC': _VISIT_SEX_AGE_BMI_CONTRA_TIME_ALCO,
        'PBMC_PER_ML': _VISIT_SEX_AGE_BMI_CONTRA_TIME_ALCO,
        'WB_PER_ML': _VISIT_SEX_AGE_BMI_CONTRA_TIME_ALCO,
        'do_not_correct': ['SAMPLE:VISIT_V2', 'SAMPLE:VISIT_V3'],
        'method': 'PCA',
        'n_components': 100,
        'visits': None,
        'get_data': get_batch_corrected_counts
    },
    'UMAP100_corrected_combat': {
        'ATAC': BLOOD + _VISIT_SEX_AGE_BMI_CONTRA_TIME_ALCO + ['RUN:TSS_ENRICHMENT'],
        'CM': WHOLE_BLOOD + _VISIT_SEX_AGE_BMI_CONTRA_TIME_ALCO,
        'CYTO': BLOOD + _VISIT_SEX_AGE_BMI_CONTRA_TIME_ALCO,
        'PBMC_PERC': _VISIT_SEX_AGE_BMI_CONTRA_TIME_ALCO,
        'PBMC_PER_ML': _VISIT_SEX_AGE_BMI_CONTRA_TIME_ALCO,
        'WB_PER_ML': _VISIT_SEX_AGE_BMI_CONTRA_TIME_ALCO,
        'do_not_correct': ['SAMPLE:VISIT_V2', 'SAMPLE:VISIT_V3'],
        'method': 'UMAP',
        'n_components': 100,
        'visits': None,
        'get_data': get_batch_corrected_counts
    },
    'PCA1000_corrected_combat': {
        'ATAC': BLOOD + _VISIT_SEX_AGE_BMI_CONTRA_TIME_ALCO + ['RUN:TSS_ENRICHMENT'],
        'CM': WHOLE_BLOOD + _VISIT_SEX_AGE_BMI_CONTRA_TIME_ALCO,
        'CYTO': BLOOD + _VISIT_SEX_AGE_BMI_CONTRA_TIME_ALCO,
        'PBMC_PERC': _VISIT_SEX_AGE_BMI_CONTRA_TIME_ALCO,
        'PBMC_PER_ML': _VISIT_SEX_AGE_BMI_CONTRA_TIME_ALCO,
        'WB_PER_ML': _VISIT_SEX_AGE_BMI_CONTRA_TIME_ALCO,
        'do_not_correct': ['SAMPLE:VISIT_V2', 'SAMPLE:VISIT_V3'],
        'method': 'PCA',
        'n_components': 1000,
        'visits': None,
        'get_data': get_batch_corrected_counts
    },
    'UMAP1000_corrected_combat': {
        'ATAC': BLOOD + _VISIT_SEX_AGE_BMI_CONTRA_TIME_ALCO + ['RUN:TSS_ENRICHMENT'],
        'CM': WHOLE_BLOOD + _VISIT_SEX_AGE_BMI_CONTRA_TIME_ALCO,
        'CYTO': BLOOD + _VISIT_SEX_AGE_BMI_CONTRA_TIME_ALCO,
        'PBMC_PERC': _VISIT_SEX_AGE_BMI_CONTRA_TIME_ALCO,
        'PBMC_PER_ML': _VISIT_SEX_AGE_BMI_CONTRA_TIME_ALCO,
        'WB_PER_ML': _VISIT_SEX_AGE_BMI_CONTRA_TIME_ALCO,
        'do_not_correct': ['SAMPLE:VISIT_V2', 'SAMPLE:VISIT_V3'],
        'method': 'UMAP',
        'n_components': 1000,
        'visits': None,
        'get_data': get_batch_corrected_counts
    },
    'combat_SX': {
        'ATAC': None, #BLOOD + _VISIT_SEX_AGE_BMI_CONTRA_TIME_ALCO + ['RUN:TSS_ENRICHMENT'],
        'CM': None,
        'CYTO': None,
        'PBMC_PERC': None,
        'PBMC_PER_ML': None,
        'WB_PER_ML': None,
        'do_not_correct': ['SAMPLE:VISIT_V2', 'SAMPLE:VISIT_V3', 'DONOR:SEX_M', 'DONOR:AGE', 'DONOR:BMI', 'DONOR:oralContraceptivesIncludingMen', 'SAMPLE:VISIT_TIME_REAL'],
        'method': None,
        'n_components': None,
        'visits': None,
        'get_data': get_batch_corrected_counts
    },
    'combat_adjScores': {
        'ATAC': None, #BLOOD + _VISIT_SEX_AGE_BMI_CONTRA_TIME_ALCO + ['RUN:TSS_ENRICHMENT'],
        'CM': None,
        'CYTO': None,
        'PBMC_PERC': None,
        'PBMC_PER_ML': None,
        'WB_PER_ML': None,
        'do_not_correct': ['SAMPLE:VISIT_V2', 'SAMPLE:VISIT_V3', 'DONOR:SEX_M', 'DONOR:AGE', 'DONOR:BMI', 'DONOR:oralContraceptivesIncludingMen', 'SAMPLE:VISIT_TIME_REAL'],
        'method': None,
        'n_components': None,
        'visits': None,
        'get_data': get_batch_corrected_counts
    },
    'combat_noS_adjScores': {
        'ATAC': None, #BLOOD + _VISIT_SEX_AGE_BMI_CONTRA_TIME_ALCO + ['RUN:TSS_ENRICHMENT'],
        'CM': None,
        'CYTO': None,
        'PBMC_PERC': None,
        'PBMC_PER_ML': None,
        'WB_PER_ML': None,
        'do_not_correct': ['SAMPLE:VISIT_V2', 'SAMPLE:VISIT_V3', 'DONOR:SEX_M', 'DONOR:AGE', 'DONOR:BMI', 'DONOR:oralContraceptivesIncludingMen', 'SAMPLE:VISIT_TIME_REAL'],
        'method': None,
        'n_components': None,
        'visits': None,
        'get_data': get_batch_corrected_counts
    },
    # 'corrected_donor': {
    #     'ATAC': BLOOD + _VISIT_DONOR_TIME + _ATAC_COVS,
    #     'CM': WHOLE_BLOOD + _VISIT_DONOR_TIME,
    #     'CYTO': BLOOD + _VISIT_DONOR_TIME,
    #     'PBMC_PERC': _VISIT_DONOR_TIME,
    #     'PBMC_PER_ML': _VISIT_DONOR_TIME,
    #     'method': None,
    #     'n_components': None,
    #     'visits': None,
    #     'get_data': get_norm_counts
    # },
    # 'corrected_combat_donor': {
    #     'ATAC': BLOOD + _VISIT_DONOR_TIME + ['RUN:TSS_ENRICHMENT'],
    #     'CM': WHOLE_BLOOD + _VISIT_DONOR_TIME,
    #     'CYTO': BLOOD + _VISIT_DONOR_TIME,
    #     'PBMC_PERC': _VISIT_DONOR_TIME,
    #     'PBMC_PER_ML': _VISIT_DONOR_TIME,
    #     'method': None,
    #     'n_components': None,
    #     'visits': None,
    #     'get_data': get_batch_corrected_counts
    # }
}

results_dir = ML_RESULTS_ROOT
print(results_dir)

model_dir = make_dir(results_dir, model)
model_fig_dir = make_dir(model_dir, 'figures')
print('\n---------\n{}\n---------\n'.format(model))

for celltype in ['PBMC']:
    print('Using this data:', MODELS[model]['get_data'])
    atac_df = MODELS[model]['get_data'](celltype).T
    annot_df = get_sample_annot(which=MODELS[model].get('which')).loc[atac_df.index]
    annot_df[['RUN:FRIP', 'RUN:ORACLE_FRIP', 'RUN:PEAKS', 'RUN:TSS_ENRICHMENT', 'RUN:FASTQC_GC_PERC']]

    if MODELS[model]['visits']:
        annot_df = annot_df.loc[annot_df['SAMPLE:VISIT'].isin(MODELS[model]['visits'])]
        atac_df = atac_df.loc[annot_df.index]
        print('Selected just these visits:', set(annot_df['SAMPLE:VISIT']))

    if MODELS[model]['ATAC']:
        print('\n--------\nCORRECTING ATAC-seq\n--------\n')
        nulls = annot_df[MODELS[model]['ATAC']].isnull().any(axis=1)
        if nulls.any():
            print(atac_df.shape)
            print('Removing {} samples with null annotations'.format(nulls.sum()))
            atac_df = atac_df.loc[~nulls]
            annot_df = annot_df.loc[~nulls]
            print(atac_df.shape)

        if model in ['corrected_LMM_combat', 'corrected_converged_LMM_combat']:
            if MODELS[model]['do_not_return_without_convergence']:
                warnings.filterwarnings('error')
            _, _, atac_df, _ = fit_linear_model(
                X_df=annot_df.loc[:, MODELS[model]['ATAC'] + ['SAMPLE:DONOR']], Y_df=atac_df,
                design=' + '.join([str(1)] + [_encode_coef(c, not is_numeric_dtype(annot_df[c].dtype)) for c in MODELS[model]['ATAC']]),
                lmm_groups=['SAMPLE:DONOR'], do_not_correct=MODELS[model]['do_not_correct'],
                return_corrected_X=True, just_correction=True, random_state=RANDOM_STATE)
            if MODELS[model]['do_not_return_without_convergence']:
                warnings.filterwarnings('default')
        else:
            atac_df, _, _ = deprecated_regress_out(
                X_df=annot_df, Y_df=atac_df,
                variables=MODELS[model]['ATAC'], do_not_correct_prefixes=MODELS[model]['do_not_correct']
            )

        print('corrected atac_df:', atac_df.shape)
        # This file is used for job_var_part_ATAC.sh
        atac_df.T.to_csv(os.path.join(DATA, 'DE', '{}_normalized_log2_CPM_PBMC.csv.gz'.format(model)))
        # sys.exit(0)

    # BLOOD HAS TO BE LAST!!!
    for d in ['CYTO', 'CM', 'PBMC_PERC', 'PBMC_PER_ML', 'WB_PER_ML']:
        if MODELS[model][d]:
            print('\n--------\nCORRECTING {}\n--------\n'.format(d))
            correct_cols = MODELS[model][d]
            assert not annot_df[correct_cols].isnull().any().any()
            _first = True
            for y in annot_df.columns[annot_df.columns.str.startswith('{}:'.format(d))]:
                nulls = annot_df[y].isnull()
                annot_df.loc[~nulls, [y]], _, _ = deprecated_regress_out(
                    X_df=annot_df.loc[~nulls], Y_df=annot_df.loc[~nulls, [y]],
                    variables=correct_cols, do_not_correct_prefixes=MODELS[model]['do_not_correct'],
                    verbose=_first
                )
                _first = False

    ax = plot_sequenced_samples(atac_df, n_samples=None, samples_as_rows=True)
    ax.set_title('{} normalized{}'.format(celltype, ' corrected' if MODELS[model]['ATAC'] is not None else ''))
    savefig(os.path.join(model_fig_dir, '{}.{}.ATAC.pdf'.format(PROJECT, celltype)))
    plt.close()

    ax = sns.scatterplot(atac_df.mean(), atac_df.std())
    ax.set_xlabel('Mean')
    ax.set_ylabel('Std')
    ax.set_title('{} normalized{}'.format(celltype, ' corrected' if MODELS[model]['ATAC'] is not None else ''))
    sns.despine()
    savefig(os.path.join(model_fig_dir, '{}.{}.ATAC_mean_std_trend.pdf'.format(PROJECT, celltype)))
    plt.close()

    # plot PCA and associations
    categorical_labels = {
            'DONOR': ['LAB:BATCH', 'SAMPLE:DONOR', 'SAMPLE:VISIT', 'DONOR:SEX', 'DONOR:AGE',
                      'SAMPLE:VISIT_TIME_REAL']
    }

    continuous_labels = {
            'PBMC_PERC': annot_df.columns[annot_df.columns.str.startswith('PBMC_PERC:')].tolist(),
            'PBMC_PER_ML': annot_df.columns[annot_df.columns.str.startswith('PBMC_PER_ML:')].tolist(),
            'WB_PER_ML': annot_df.columns[annot_df.columns.str.startswith('WB_PER_ML:')].tolist(),
            'CYTO': annot_df.columns[annot_df.columns.str.contains('^CYTO:.*_good$')].tolist(),
            'CM': annot_df.columns[annot_df.columns.str.startswith('CM:')].tolist(),
            'RUN': ['RUN:FRIP', 'RUN:ORACLE_FRIP', 'RUN:PEAKS', 'RUN:TSS_ENRICHMENT', 'RUN:FASTQC_GC_PERC']
    }

    pca_analysis(
        categorical_labels, continuous_labels,
        df=atac_df,
        annot_df=annot_df,
        title_prefix=celltype,
        fig_dir=model_fig_dir, fig_prefix='{}.{}.{}.ATAC_PCA.'.format(PROJECT, celltype, model),
        random_state=RANDOM_STATE, expl_var_thr=21,
        scatter_kwargs=dict(n_dims=5, palettes=None, legends={'LAB:BATCH': None, 'SAMPLE:DONOR': None}, orders='random',
                            num_legend_fmt='.1f', alpha=0.5, s=50),
        categorical_heatmap_kwargs=dict(n_dims=10, vmin=0, vmax=0.2, center=None, annot=True, fmt='.0e', annot_kws={'size': 10}),
        continuous_heatmap_kwargs=dict(n_dims=10, vmin=-0.5, vmax=0.5, center=0, annot=True, fmt='.2f', annot_kws={'size': 10}),
        plot_variance_explained=True, plot_scatter=False, plot_assoc_heatmap=True, verbose=False,
        show=False, save_fig=True, fig_ext='pdf'
    )

    if SAVE_DATA:
        if MODELS[model]['n_components']:
            print('Transforming atac_df with PCA n_components', MODELS[model]['n_components'])
            pca = PCA(n_components=MODELS[model]['n_components'], random_state=RANDOM_STATE)
            atac_df = pd.DataFrame(data=pca.fit_transform(StandardScaler().fit_transform(atac_df)),
                                  index=atac_df.index)
            atac_df.columns = ['PC{}'.format(i + 1) for i in range(atac_df.shape[1])]
            print(atac_df.shape)
        else:
            de_peaks_df = pd.DataFrame()
            # get differential accessibility results
            # fn = de_fn(celltype, model='donor_as_mixed.batch.sex.age.blood.TSS_enr.visit_time.thm.innate_nonspecific_24h_wo_LAC_IL10_IL1ra_V3_FC1.2_responder')
            # de_peaks_df = pd.read_csv(fn, index_col=0)
            # de_peaks_df = de_peaks_df.loc[atac_df.columns]
            # prefix = 'thm.innate_nonspecific_24h_V3_FC1.2_'
            # new_column_names = {'A': 'intensity'}
            # for c in ['N', 'N.V2', 'N.V3', 'R.V2', 'R.V3']:
            #     for stat_old, stat_new in [('Coef', 'LFC'), ('p.value', 'pval')]:
            #         new_column_names['{}.{}{}'.format(stat_old, prefix, c)] = '{}_{}'.format(stat_new, c)
            # de_peaks_df = de_peaks_df.rename(new_column_names, axis=1)[new_column_names.values()]
            # print(de_peaks_df.columns)
            de_peaks_df = pd.DataFrame()

            # get peaks annotations and put it all together
            peaks_info_df = get_peak_annot().loc[atac_df.columns]
            peak_annot_df = pd.concat([de_peaks_df, peaks_info_df], axis=1)
            assert peak_annot_df.index.equals(atac_df.columns)

        assert annot_df.index.equals(atac_df.index)
        annot_df = annot_df.astype({c: 'float' for c in NAN_BOOLS_ANNOT_COLS})

        # save
        data_fn = os.path.join(model_dir, '{}.{}.hdf'.format(PROJECT, model))
        for key, df in [('{}/atac'.format(celltype), atac_df),
                        ('{}/sample_annot'.format(celltype), annot_df)] + \
                       ([] if MODELS[model]['n_components'] else [('{}/peak_annot'.format(celltype), peak_annot_df)]):
            with warnings.catch_warnings():
                warnings.simplefilter(action='ignore', category=tables.NaturalNameWarning)
                warnings.simplefilter(action='ignore', category=pd.errors.PerformanceWarning)
                df.to_hdf(data_fn, key=key, complevel=9)
