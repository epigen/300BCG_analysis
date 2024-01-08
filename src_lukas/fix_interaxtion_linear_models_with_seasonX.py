import matplotlib as mpl
mpl.use('Agg')
import misc
import bcg_utils as utils
import os
import pandas as pd
import numpy as np
from linear_models import _Y_blood_types, Y_regex, Y_transform, BASELINE_COLS


def main():
    celltype = 'PBMC'
    for Y_name, IC_blood in [
        # ('CYTO', True),
        ('CYTO', False),
        # ('CM', False),
        # ('WB_PER_ML', False),
        # ('WB_PERC', False)
    ]:
        print('\n', Y_name, '\n')
        if IC_blood:
            COEFS = [(['{}:{}'.format('WB_PER_ML' if Y_name == 'CYTO' else 'WB_PER_ML' if Y_name == 'CM' else 'XX', b)],
                      ['IC_V1_CORR_WB_PER_ML:{}'.format(b)]) for b in _Y_blood_types]
        else:
            COEFS = [
                # (['SAMPLE:VISIT_DATE_2PI_SIN', 'SAMPLE:VISIT_DATE_2PI_COS'], ['DONOR:IC_DATE_2PI_SIN', 'DONOR:IC_DATE_2PI_COS']),
                (['SAMPLE:VISIT_TIME_REAL'], ['DONOR:IC_TIME_REAL']),
                (['SAMPLE:alcoholInLast24h'], ['DONOR:IC_alcoholInLast24h'])
            ]

        results_df = []
        for sample_coefs, V1_IC_coefs in COEFS:
            annot_df = misc.get_sample_annot()
            annot_df = annot_df.rename({'CYTO:C.albicans.yeast_24h_PBMC_IFNg_excluded': 'CYTO:C.albicans.yeast_24h_PBMC_IFNg_good'}, axis=1)
            annot_df = annot_df.loc[(annot_df['SAMPLE:TISSUE'] == celltype) & (~annot_df['DONOR:IC_EVENING'])]
            Y_df = annot_df.loc[:, annot_df.columns.str.contains(Y_regex[Y_name])].copy()
            print(Y_name, Y_df.shape)
            if Y_transform[Y_name]:
                print('Transforming', Y_name, Y_transform[Y_name])
                Y_df = Y_transform[Y_name](Y_df)
            cols = sorted(set(['SAMPLE:DONOR', 'SAMPLE:VISIT', 'DONOR:IC_DATE_2PI_SIN', 'DONOR:IC_DATE_2PI_COS'] + BASELINE_COLS[Y_name] + sample_coefs + V1_IC_coefs \
                              + (misc.BLOOD if Y_name == 'CYTO' and not IC_blood else misc.WHOLE_BLOOD if Y_name == 'CM' and not IC_blood else [])))
            X_df = annot_df[cols].copy()
            X_df, Y_df = misc.drop_na_for_LM(X_df, Y_df)
            X_df = pd.concat([pd.get_dummies(X_df[[c for c in X_df.columns if c != 'SAMPLE:DONOR']], drop_first=True),
                              X_df['SAMPLE:DONOR']], axis=1)
            X_df.columns = utils.replace_chars_for_LM(X_df.columns)
            for v in ['V2', 'V3']:
                for c in ['DONOR_IC_DATE_2PI_SIN', 'DONOR_IC_DATE_2PI_COS']:
                    X_df['{}_{}'.format(v, c)] = X_df['SAMPLE_VISIT_{}'.format(v)] * X_df[c]
            print(X_df.columns)
            print(X_df)

            contrasts, do_not_correct = [], ['SAMPLE_VISIT_V2', 'SAMPLE_VISIT_V3']
            for v in ['V2', 'V3']:
                _con = []
                for c in misc.replace_chars_for_LM(pd.Series(V1_IC_coefs)):
                    _interaction = '{}_{}'.format(v, c)
                    _con.append(_interaction)
                    do_not_correct.append(_interaction)
                    X_df[_interaction] = X_df['SAMPLE_VISIT_{}'.format(v)] * X_df[c]
                contrasts.append(_con)
            X_df = X_df.where(X_df != 0, other=0)
            print('contrasts:', contrasts)
            print('do_not_correct:', do_not_correct)

            # DROP THE EXTRA NON-SENSICAL CORRECTION
            X_df = X_df.drop(misc.replace_chars_for_LM(pd.Series(V1_IC_coefs)), axis=1)
            design = '1 + ' + misc.replace_chars_for_LM(' + '.join(
                [utils._encode_coef(c, False) for c in X_df.columns if c != 'SAMPLE_DONOR']))
            print('design:', design)

            _, coefs_df, contrasts_df, corrected_df, _ = utils.fit_linear_model(
                X_df, Y_df, design=design, lmm_groups=['SAMPLE_DONOR'], contrasts=contrasts,
                do_not_correct=do_not_correct, return_corrected_X=True,
                random_state=misc.RANDOM_STATE, verbose=True)

            for c in contrasts:
                if len(c) == 1:
                    c_mask = coefs_df.index.get_level_values('contrast') == c[0]
                    coefs_df.loc[c_mask, 'padj'] = misc.adjusted_pvals(coefs_df.loc[c_mask, 'p.value'])
                    results_df.append(coefs_df.loc[c_mask])
                else:
                    c_mask = contrasts_df.index.get_level_values('contrast') == '__'.join(c)
                    contrasts_df.loc[c_mask, 'padj'] = misc.adjusted_pvals(contrasts_df.loc[c_mask, 'p.value'])
                    contrasts_df.loc[c_mask, 'Coef'] = np.nan
                    results_df.append(contrasts_df.loc[c_mask])

            corrected_df.to_csv(
                    os.path.join('results', 'LR', '{}.seasonX.{}.{}'.format(celltype, Y_name, (
                        'BLOOD' if Y_name == 'CYTO' and not IC_blood else 'WHOLE_BLOOD' if Y_name == 'CM' and not IC_blood else 'None')),
                                 '{}.corrected.fixed_{}_interactions.{}.csv'.format(Y_name, 'factor' if not IC_blood else 'blood', '__'.join(V1_IC_coefs).replace(os.path.sep, '_'))))

            if len(V1_IC_coefs) > 1:
                coefs_df.index = pd.MultiIndex.from_arrays(
                    [coefs_df.index.get_level_values('target'),
                     coefs_df.index.get_level_values('contrast').str.replace(
                         '^V2_', 'SAMPLE_VISIT[T.V2]:').str.replace('^V3_', 'SAMPLE_VISIT[T.V3]:')])
                coefs_df.to_csv(
                    os.path.join('results', 'LR', '{}.seasonX.{}.{}'.format(celltype, Y_name, (
                        'BLOOD' if Y_name == 'CYTO' and not IC_blood else 'WHOLE_BLOOD' if Y_name == 'CM' and not IC_blood else 'None')),
                                 '{}.coefs.fixed_{}_interactions.{}.csv'.format(Y_name, 'factor' if not IC_blood else 'blood', '__'.join(V1_IC_coefs).replace(os.path.sep, '_'))))

        results_df = pd.concat(results_df)
        results_df.index = pd.MultiIndex.from_arrays(
            [results_df.index.get_level_values('target'),
             results_df.index.get_level_values('contrast').str.replace(
                 '^V2_', 'SAMPLE_VISIT[T.V2]:').str.replace('^V3_', 'SAMPLE_VISIT[T.V3]:').str.replace(
                 '__V2_', '__SAMPLE_VISIT[T.V2]:').str.replace('__V3_', '__SAMPLE_VISIT[T.V3]:')]
        )
        results_df[['Coef', 'p.value', 'padj', 'resid.df', 'stat']].to_csv(
            os.path.join('results', 'LR', '{}.seasonX.{}.{}'.format(celltype, Y_name, ('BLOOD' if Y_name == 'CYTO' and not IC_blood else 'WHOLE_BLOOD' if Y_name == 'CM' and not IC_blood else 'None')),
                         '{}.fixed_{}_interactions.csv'.format(Y_name, 'factor' if not IC_blood else 'blood')))
        print('\n----------------------------\n')

    # sample_coefs = ['WB_PER_ML:T/CD4/TREG']
    # V1_IC_coefs = ['IC_WB_PER_ML:T/CD4/TREG']
    #
    # visit = 'V1'
    # annot_df = misc.get_sample_annot()
    # annot_df = annot_df.loc[(annot_df['SAMPLE:TISSUE'] == celltype) & (annot_df['SAMPLE:VISIT'] == visit)]
    # cyto_df = annot_df.loc[:, annot_df.columns.str.contains('^CYTO:.*_good$')].copy()
    # cols = sorted(set(['DONOR:AGE', 'DONOR:SEX', 'SAMPLE:VISIT_TIME_REAL'] + misc.BLOOD + sample_coefs))
    # X_df = annot_df[cols].astype({c: float for c in sample_coefs}).copy()
    # X_df, cyto_df = misc.drop_na_for_LM(X_df, cyto_df)
    # _, V1_coefs_df, V1_contrasts_df, V1_corr_cyto_df, _ = utils.fit_linear_model(
    #     X_df, cyto_df,
    #     design='1 + ' + misc.replace_chars_for_LM(' + '.join(X_df.columns)),
    #     contrasts=[misc.replace_chars_for_LM(pd.Series(sample_coefs)).tolist()],
    #     do_not_correct=misc.replace_chars_for_LM(pd.Series(sample_coefs)).tolist(),
    #     return_corrected_X=True, random_state=misc.RANDOM_STATE, verbose=True)
    # V1_coefs_df['padj'] = misc.adjusted_pvals(V1_coefs_df['p.value'])
    #
    # for visit_interaction in ['V3']:
    #     _sample_corrections = sorted(set(['SAMPLE:VISIT_TIME_REAL'] + misc.BLOOD + sample_coefs))
    #     annot_df = misc.get_sample_annot()
    #     annot_df = annot_df.loc[(annot_df['SAMPLE:TISSUE'] == celltype) & (annot_df['DONOR:IC_EVENING'] == False)]
    #     cyto_df = annot_df.loc[:, annot_df.columns.str.contains('^CYTO:.*_good$')].copy()
    #     X_df = annot_df[_sample_corrections].astype({c: float for c in sample_coefs}).copy()
    #     X_df, cyto_df = misc.drop_na_for_LM(X_df, cyto_df)
    #     cyto_df = cyto_df - np.dot(X_df[_sample_corrections], V1_coefs_df['Coef'].unstack()[
    #         misc.replace_chars_for_LM(pd.Series(_sample_corrections))].T)
    #
    #     cyto_df.index = cyto_df.index.str.split('_', expand=True)
    #     cyto_df.index.names = ['SAMPLE:DONOR', 'SAMPLE:VISIT', 'SAMPLE:TISSUE']
    #     cyto_df = misc.visit_diff(cyto_df, visit=visit_interaction, base_visit='V1')
    #     cyto_df.index = cyto_df.index + '_{}_PBMC'.format(visit_interaction)
    #     cyto_df.index.name = annot_df.index.name
    #     cyto_df = cyto_df.loc[~cyto_df.isnull().all(axis=1)]
    #     annot_df = annot_df.loc[cyto_df.index].sort_values(['DONOR:IC_DATE_REAL', 'SAMPLE:VISIT', 'SAMPLE:DONOR'])
    #     cyto_df = cyto_df.loc[annot_df.index]
    #
    #     X_df = annot_df[V1_IC_coefs].astype({c: float for c in V1_IC_coefs}).copy()
    #     X_df, cyto_df = misc.drop_na_for_LM(X_df, cyto_df)
    #     _, coefs_df, contrasts_df, corr_cyto_df, _ = utils.fit_linear_model(
    #         X_df, cyto_df,
    #         design='1 + ' + misc.replace_chars_for_LM(' + '.join(X_df.columns)),
    #         contrasts=[misc.replace_chars_for_LM(pd.Series(V1_IC_coefs)).tolist()],
    #         do_not_correct=misc.replace_chars_for_LM(pd.Series(V1_IC_coefs)).tolist(),
    #         return_corrected_X=True, random_state=misc.RANDOM_STATE, verbose=True)
    #     contrasts_df['padj'] = misc.adjusted_pvals(contrasts_df['p.value'])
    #     print(visit_interaction)
    #     print(contrasts_df.sort_values('p.value')['p.value'].droplevel('contrast').head(10))


if __name__ == '__main__':
    main()
