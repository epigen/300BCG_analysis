import sys
import misc
import bcg_utils as utils
import os
import pandas as pd
import numpy as np
from linear_models import _Y_blood_types, Y_regex, Y_transform, BASELINE_COLS

results_dir = 'results'
IC_BLOOD = False
CELLTYPE = 'PBMC'


def main():
    Y_name = sys.argv[1]
    assert Y_name in ['CYTO', 'CM', 'WB_PER_ML']
    print('\n', Y_name, '\n')

    if IC_BLOOD:
        COEFS = [(['{}:{}'.format('WB_PER_ML' if Y_name == 'CYTO' else 'WB_PER_ML' if Y_name == 'CM' else 'XX', b)],
                  ['IC_V1_CORR_WB_PER_ML:{}'.format(b)]) for b in _Y_blood_types]
    else:
        COEFS = [
            (['SAMPLE:VISIT_DATE_2PI_SIN', 'SAMPLE:VISIT_DATE_2PI_COS'], ['DONOR:IC_DATE_2PI_SIN', 'DONOR:IC_DATE_2PI_COS']),
            (['SAMPLE:VISIT_TIME_REAL'], ['DONOR:IC_TIME_REAL']),
            (['SAMPLE:alcoholInLast24h'], ['DONOR:IC_alcoholInLast24h'])
        ]

    for sample_coefs, V1_IC_coefs in COEFS:
        annot_df = misc.get_sample_annot()
        annot_df = annot_df.loc[(annot_df['SAMPLE:TISSUE'] == CELLTYPE) & (~annot_df['DONOR:IC_EVENING'])]
        Y_df = annot_df.loc[:, annot_df.columns.str.contains(Y_regex[Y_name])].copy()
        print(Y_name, Y_df.shape)
        if Y_transform[Y_name]:
            print('Transforming', Y_name, Y_transform[Y_name])
            Y_df = Y_transform[Y_name](Y_df)
        cols = sorted(set(['SAMPLE:DONOR', 'SAMPLE:VISIT'] + BASELINE_COLS[Y_name] + sample_coefs + V1_IC_coefs \
                          + (misc.BLOOD if Y_name == 'CYTO' and not IC_BLOOD else misc.WHOLE_BLOOD if Y_name == 'CM' and not IC_BLOOD else [])))
        X_df = annot_df[cols].copy()
        X_df, Y_df = misc.drop_na_for_LM(X_df, Y_df)
        X_df = pd.concat([pd.get_dummies(X_df[[c for c in X_df.columns if c != 'SAMPLE:DONOR']], drop_first=True),
                          X_df['SAMPLE:DONOR']], axis=1)
        X_df.columns = utils.replace_chars_for_LM(X_df.columns)

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

        results_df = {}
        todo = []
        for c in contrasts:
            if len(c) == 1:
                c_mask = coefs_df.index.get_level_values('contrast') == c[0]
                coefs_df.loc[c_mask, 'padj'] = misc.adjusted_pvals(coefs_df.loc[c_mask, 'p.value'])
                results_df['__'.join(c)] = coefs_df.loc[c_mask]
                results_df['__'.join(c)].index = pd.MultiIndex.from_arrays(
                    [results_df['__'.join(c)].index.get_level_values('target'),
                     results_df['__'.join(c)].index.get_level_values('contrast').str.replace(
                         '^V2_', 'SAMPLE_VISIT[T.V2]:').str.replace('^V3_', 'SAMPLE_VISIT[T.V3]:').str.replace(
                         '__V2_', '__SAMPLE_VISIT[T.V2]:').str.replace('__V3_', '__SAMPLE_VISIT[T.V3]:')])
                res_fn = misc.summary_lr_fn(celltype=CELLTYPE, model=f'final.{Y_name}.{"BLOOD" if Y_name == "CYTO" else "WHOLE_BLOOD" if Y_name == "CM" else "None"}',
                                            Y_name=Y_name, visit_interaction=True, LMM=True, remove_evening=True, results_dir=results_dir)
                res_df = pd.read_csv(res_fn, index_col=['target', 'contrast'])
                res_df.loc[results_df['__'.join(c)].index] = results_df['__'.join(c)].loc[:, res_df.columns]
                res_df.to_csv(res_fn)
            else:
                todo.append(c)
                c_mask = contrasts_df.index.get_level_values('contrast') == '__'.join(c)
                contrasts_df.loc[c_mask, 'padj'] = misc.adjusted_pvals(contrasts_df.loc[c_mask, 'p.value'])
                contrasts_df.loc[c_mask, 'Coef'] = np.nan
                results_df['__'.join(c)] = contrasts_df.loc[c_mask]
                results_df['__'.join(c)].index = pd.MultiIndex.from_arrays(
                    [results_df['__'.join(c)].index.get_level_values('target'),
                     results_df['__'.join(c)].index.get_level_values('contrast').str.replace(
                         '^V2_', 'SAMPLE_VISIT[T.V2]:').str.replace('^V3_', 'SAMPLE_VISIT[T.V3]:').str.replace(
                         '__V2_', '__SAMPLE_VISIT[T.V2]:').str.replace('__V3_', '__SAMPLE_VISIT[T.V3]:')]
                )

        if len(todo) != 0:
            formatted_coefs = f'.{"__".join(V1_IC_coefs).replace(os.path.sep, "_")}'
            fn_template = os.path.join(
                results_dir, 'LR',
                f'{CELLTYPE}.final.{Y_name}.{"BLOOD" if Y_name == "CYTO" and not IC_BLOOD else "WHOLE_BLOOD" if Y_name == "CM" and not IC_BLOOD else "None"}',
                f'LR_{"{data}"}_{Y_name}{"{prefix}"}.visit_interaction.LMM.remove_evening{".IC_blood" if IC_BLOOD else ""}{"{suffix}"}.csv')

            pd.concat([results_df['__'.join(c)] for c in todo]).to_csv(
                fn_template.format(data='results', prefix='', suffix=formatted_coefs))

            coefs_df.index = pd.MultiIndex.from_arrays(
                [coefs_df.index.get_level_values('target'),
                 coefs_df.index.get_level_values('contrast').str.replace(
                     '^V2_', 'SAMPLE_VISIT[T.V2]:').str.replace('^V3_', 'SAMPLE_VISIT[T.V3]:').str.replace(
                     '__V2_', '__SAMPLE_VISIT[T.V2]:').str.replace('__V3_', '__SAMPLE_VISIT[T.V3]:')])
            coefs_df.to_csv(fn_template.format(data='coefs', prefix=formatted_coefs, suffix=''))

            corrected_df.to_csv(fn_template.format(data='corrected', prefix=formatted_coefs, suffix=''))


if __name__ == '__main__':
    main()
