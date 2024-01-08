import matplotlib as mpl
mpl.use('Agg')
import argparse
import misc
import bcg_utils as utils
import os
import sys
import traceback
import pandas as pd
import numpy as np
from pandas.api.types import is_numeric_dtype


def _is_single_visit(visits):
    return visits is not None and (isinstance(visits, str) or len(visits) == 1)


def _correct_Y(X_df, Y_df, design, lmm_groups, do_not_correct):
    _X_df, _Y_df = misc.drop_na_for_LM(X_df.loc[Y_df.index], Y_df, verbose=True)
    _, _, _Y_df, _ = utils.fit_linear_model(
        _X_df, _Y_df, design=design, lmm_groups=lmm_groups, do_not_correct=do_not_correct, return_corrected_X=True,
        random_state=misc.RANDOM_STATE, verbose=True)
    Y_df.loc[:, :] = np.nan
    Y_df.loc[_Y_df.index, Y_df.columns] = _Y_df.loc[:, Y_df.columns]
    return Y_df


_Y_blood_types = ['MONO/CLASSICAL', 'MONO/INTERMEDIATE', 'MONO/NON_CLASSICAL', 'T/CD8', 'T/CD4', 'T/CD4/TREG', 'B', 'NK', 'NKT', 'NEUTRO']
Y_regex = {
    'CYTO': '^CYTO:.*_good$',
    'CM': '^CM:',
    'WB_PER_ML': '|'.join(['^WB_PER_ML:{}$'.format(b) for b in _Y_blood_types]),
}

Y_transform = {
    'CYTO': None,
    'CM': None,
    'WB_PER_ML': np.log2,
}

_age_sex_time = ['DONOR:SEX', 'DONOR:AGE', 'SAMPLE:VISIT_TIME_REAL']
BASELINE_COLS = {
    'CYTO': _age_sex_time,
    'CM': _age_sex_time + ['DONOR:BMI', 'DONOR:oralContraceptivesIncludingMen'],
    'WB_PER_ML': _age_sex_time,
}

BLOOD_NAMES = {
    'None': [],
    'BLOOD': misc.BLOOD,
    'WHOLE_BLOOD': misc.WHOLE_BLOOD,
}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--config')
    parser.add_argument('--celltype', default='PBMC')
    parser.add_argument('--model', default='test')
    parser.add_argument('--annot_fn', default=misc.get_sample_annot_fn())
    parser.add_argument('--Y_name', default='CYTO', choices=Y_regex.keys())
    parser.add_argument('--blood', choices=BLOOD_NAMES, default="None")
    parser.add_argument('--extras', type=str, nargs='*')
    parser.add_argument('--coefs', type=str, nargs='*')
    parser.add_argument('--include_host_factors', action='store_true')
    parser.add_argument('--include_season', action='store_true')
    parser.add_argument('--include_blood_perc', action='store_true')
    parser.add_argument('--include_blood_per_ml', action='store_true')
    parser.add_argument('--include_CM', action='store_true')
    parser.add_argument('--include_cyto', action='store_true')
    parser.add_argument('--include_scores', action='store_true')
    parser.add_argument('--visits', type=str, nargs='*')
    parser.add_argument('--V1', action='store_true')
    parser.add_argument('--V2', action='store_true')
    parser.add_argument('--V3', action='store_true')
    parser.add_argument('--visit_interaction', action='store_true') #, choices=['V2', 'V3'])
    parser.add_argument('--visit_interaction_and_correct_evening', action='store_true') #, choices=['V2', 'V3'])
    parser.add_argument('--visit_interaction_and_correct_season', action='store_true') #, choices=['V2', 'V3'])
    parser.add_argument('--LMM', action='store_true')
    parser.add_argument('--scale', action='store_true')
    parser.add_argument('--steady_state_season_hack', action='store_true')
    parser.add_argument('--donor_fixed_effects', action='store_true')
    parser.add_argument('--remove_evening', action='store_true')
    parser.add_argument('--correct_evening', action='store_true')
    parser.add_argument('--correct_season', action='store_true')
    parser.add_argument('--fold_changes', choices=['V2', 'V3'])
    parser.add_argument('--V2_fold_changes', action='store_true')
    parser.add_argument('--V3_fold_changes', action='store_true')
    parser.add_argument('--fdr', type=float, default=0.05)
    parser.add_argument('--results', default='results')
    parser.add_argument('--save_outputs', action='store_true')
    parser.add_argument('--save_models', action='store_true')
    parser.add_argument('--return_corrected_X', action='store_true')
    parser.add_argument('--correct_blood_online', action='store_true')
    parser.add_argument('--variance_partition', action='store_true')
    parser.add_argument('--verbose', type=int, default=1)

    # parser.add_argument('--figsize', type=float, nargs='+', default=[5, 4.5])
    # parser.add_argument('--dpi', type=int, default=300)
    # parser.add_argument('--not_rasterized', action='store_true')
    # parser.add_argument('--plot_context', default='talk', choices=utils.PLOT_CONTEXTS.keys())
    # parser.add_argument('--font_scale', type=float, default=1)
    # parser.add_argument('--fig_style', default='ticks')

    args = parser.parse_args()
    assert np.sum([args.V1, args.V2, args.V3]) < 2
    if args.V1:
        assert args.visits is None
        args.visits = ['V1']
    elif args.V2:
        assert args.visits is None
        args.visits = ['V2']
    elif args.V3:
        assert args.visits is None
        args.visits = ['V3']

    assert np.sum([args.V2_fold_changes, args.V3_fold_changes]) < 2
    if args.V2_fold_changes:
        assert args.fold_changes is None
        args.fold_changes = 'V2'
    elif args.V3_fold_changes:
        assert args.fold_changes is None
        args.fold_changes = 'V3'
    # assert not args.figsize or len(args.figsize) == 2
    assert not (args.config and (args.model or args.celltype))

    assert not (args.fold_changes and args.LMM)
    assert not (args.fold_changes and args.visit_interaction)
    assert not (args.fold_changes and args.visit_interaction_and_correct_evening)
    assert not (args.visit_interaction and args.visit_interaction_and_correct_evening)
    assert not (args.fold_changes and args.visit_interaction_and_correct_season)
    assert not (args.visit_interaction and args.visit_interaction_and_correct_season)
    assert not (args.visit_interaction_and_correct_evening and args.visit_interaction_and_correct_season)

    if 'SLURM_ARRAY_TASK_ID' in os.environ:
        permutation = int(os.environ['SLURM_ARRAY_TASK_ID'])
        raise ValueError
    else:
        permutation = None

    utils.make_dir(args.results, 'LR', '{}.{}'.format(args.celltype, args.model))
    fn_suffix = misc.lr_fn_suffix(visits=args.visits, visit_interaction=args.visit_interaction,
                                  visit_interaction_and_correct_evening=args.visit_interaction_and_correct_evening, visit_interaction_and_correct_season=args.visit_interaction_and_correct_season, fold_changes=args.fold_changes,
                                  LMM=args.LMM, scale=args.scale, remove_evening=args.remove_evening, correct_evening=args.correct_evening, correct_season=args.correct_season)
    fn = misc.summary_lr_fn(celltype=args.celltype, model=args.model, Y_name=args.Y_name, data='{data}',
                            visits=args.visits, visit_interaction=args.visit_interaction,
                            visit_interaction_and_correct_evening=args.visit_interaction_and_correct_evening,
                            visit_interaction_and_correct_season=args.visit_interaction_and_correct_season,
                            fold_changes=args.fold_changes, LMM=args.LMM, scale=args.scale,
                            include_season=args.include_season,
                            include_cyto=args.include_cyto, include_CM=args.include_CM,
                            include_blood_perc=args.include_blood_perc, include_blood_per_ml=args.include_blood_per_ml,
                            include_scores=args.include_scores,
                            remove_evening=args.remove_evening,
                            correct_evening=args.correct_evening,
                            correct_season=args.correct_season, extras='__'.join(args.coefs) if args.coefs else None,
                            results_dir=args.results)
    if os.path.exists(fn.format(data='results')):
        print('Already finished.', file=sys.stderr)
    else:
        print('Running analyses for', fn.format(data='results'), file=sys.stderr)
        results_df, full_results_df = run_all_coefs(celltype=args.celltype, model=args.model, Y_name=args.Y_name, coefs=args.coefs,
                                   include_season=args.include_season, include_host_factors=args.include_host_factors, include_blood_perc=args.include_blood_perc, include_blood_per_ml=args.include_blood_per_ml, include_cyto=args.include_cyto, include_CM=args.include_CM,
                                   include_scores=args.include_scores,
                                visits=args.visits, blood=BLOOD_NAMES[args.blood], extras=args.extras, annot_fn=args.annot_fn,
                                visit_interaction=args.visit_interaction, visit_interaction_and_correct_evening=args.visit_interaction_and_correct_evening,
                                   visit_interaction_and_correct_season=args.visit_interaction_and_correct_season,
                                   LMM=args.LMM, scale=args.scale, steady_state_season_hack=args.steady_state_season_hack,
                                   donor_fixed_effects=args.donor_fixed_effects, remove_evening=args.remove_evening,
                                   variance_partition=args.variance_partition, correct_blood_online=args.correct_blood_online,
                                   correct_evening=args.correct_evening, correct_season=args.correct_season,
                                fold_changes=('SAMPLE:VISIT', 'V1', args.fold_changes) if args.fold_changes else None,
                                permutation=permutation,
                                fdr=args.fdr,
                                return_corrected_X=args.return_corrected_X, save_outputs=args.save_outputs, save_models=args.save_models, fn_suffix=fn_suffix, verbose=args.verbose)
        results_df.to_csv(fn.format(data='results'))
        full_results_df.to_csv(fn.format(data='full_results'))
        print('Done!')


def run_all_coefs(celltype, model, Y_name, coefs, include_season, include_host_factors, include_blood_perc, include_blood_per_ml, include_cyto, include_CM, include_scores, visits, blood, extras, annot_fn, visit_interaction, visit_interaction_and_correct_evening, visit_interaction_and_correct_season,
                  LMM, scale, steady_state_season_hack, donor_fixed_effects, remove_evening, variance_partition, correct_blood_online, correct_evening, correct_season, fold_changes, permutation, fdr, return_corrected_X, save_outputs, save_models, fn_suffix, verbose):
    if verbose:
        print(annot_fn)
    annot_df = misc.get_sample_annot(fn=annot_fn)
    assert not (fold_changes and LMM)
    assert not (fold_changes and visit_interaction)
    assert not (fold_changes and visit_interaction_and_correct_evening)
    assert not (fold_changes and visit_interaction_and_correct_season)

    fold_change_or_interatcion_effect = visit_interaction or visit_interaction_and_correct_evening or visit_interaction_and_correct_season or fold_changes

    if coefs:
        COEFS = [coefs]
    elif include_host_factors:
        COEFS = pd.read_csv(os.path.join(misc.METADATA, 'selected_factors.csv'), comment='#', header=None)[0]
        if fold_change_or_interatcion_effect:
            COEFS = COEFS[COEFS.str.contains('^DONOR:')].tolist()
        else:
            COEFS = COEFS[~COEFS.str.contains('^DONOR:IC_')].tolist()
    elif include_season:
        COEFS = ([['DONOR:IC_DATE_2PI_SIN', 'DONOR:IC_DATE_2PI_COS']] if fold_change_or_interatcion_effect else [['SAMPLE:VISIT_DATE_2PI_SIN', 'SAMPLE:VISIT_DATE_2PI_COS']])
    elif include_scores:
        COEFS = annot_df.columns[annot_df.columns.str.contains('^thm.adaptive_MTB_7d_V3|^thm.heterologous_nonspecific_7d_V3|^thm.innate_nonspecific_24h_wo_LAC_IL10_IL1ra_V3|^DONOR:scarSize_v3|^DONOR:sizeOfVaccinationSpot_v2')].tolist()
    elif include_blood_perc:
        if fold_change_or_interatcion_effect:
            COEFS = ['{}:{}'.format('IC_CORR_PBMC_PERC' if Y_name.startswith('CORR_') else 'IC_PBMC_PERC', ':'.join(b.split(':')[1:])) for b in BLOOD_NAMES['BLOOD']]
        else:
            COEFS = list(BLOOD_NAMES['BLOOD'])
    elif include_blood_per_ml:
        if fold_change_or_interatcion_effect:
            COEFS = ['{}:{}'.format('IC_CORR_WB_PER_ML' if Y_name.startswith('CORR_') else 'IC_WB_PER_ML', b) for b in _Y_blood_types]
        else:
            COEFS = ['WB_PER_ML:{}'.format(b) for b in _Y_blood_types]
    elif include_CM:
        COEFS = annot_df.columns[annot_df.columns.str.contains('^IC_CM:')].tolist()
    elif include_cyto:
        COEFS = annot_df.columns[annot_df.columns.str.contains('^IC_CYTO:.*_good')].tolist()
    else:
        COEFS = None
    assert COEFS is None or np.isin([c for coefs in COEFS for c in (coefs if utils.is_iterable(coefs) else [coefs])], annot_df.columns).all(), COEFS
    if verbose:
        print(COEFS, '\n')

    results, full_results = [], []
    for new_coefs in (COEFS if COEFS is not None else [COEFS]):
        if isinstance(new_coefs, str):
            new_coefs = [new_coefs]
        try:
            if verbose:
                print(new_coefs)
            if not (_is_single_visit(visits) and fold_change_or_interatcion_effect):
                results_df, full_results_df = run_linear_models(celltype=celltype, model=model, Y_name=Y_name, new_coefs=new_coefs,
                                              new_coef_dtype=annot_df[new_coefs[0]].dtype if new_coefs is not None and len(new_coefs) == 1 else None,
                                              visits=visits,
                                              blood=blood, extras=(extras if extras else []) + \
                                                                  (['DONOR:IC_EVENING'] if correct_evening else []) + \
                                                                  (['DONOR:IC_DATE_2PI_SIN', 'DONOR:IC_DATE_2PI_COS'] if correct_season else []),
                                              annot_fn=annot_fn, annot_df=annot_df,
                                              visit_interaction=visit_interaction, visit_interaction_and_correct_evening=visit_interaction_and_correct_evening, visit_interaction_and_correct_season=visit_interaction_and_correct_season, LMM=LMM, scale=scale, donor_fixed_effects=donor_fixed_effects,
                                              variance_partition=variance_partition, correct_blood_online=correct_blood_online,
                                              remove_samples=annot_df.loc[annot_df['DONOR:IC_EVENING']].index.tolist() if remove_evening else None,
                                              fold_changes=fold_changes, steady_state_season_hack=steady_state_season_hack, permutation=permutation, fdr=fdr,
                                              return_corrected_X=return_corrected_X, save_outputs=save_outputs, save_models=save_models, fn_suffix='{}.{}'.format('__'.join(new_coefs) if new_coefs is not None else None, fn_suffix), verbose=verbose)

                results.append(results_df)
                full_results.append(full_results_df)
                if verbose and 'padj' in results_df.columns:
                    print((results_df['padj'] < fdr).astype(int).groupby('contrast').sum())
            else:
                print('Skipping because interaction and fold-changes cannot be done with a single visit\n')
        except Exception:
            print('\nERROR with', new_coefs)
            traceback.print_exception(*sys.exc_info())
            print('\n')
        sys.stdout.flush()

    results_df = pd.concat(results, axis=0).sort_index(axis=1)
    full_results_df = pd.concat(full_results, axis=0).sort_index(axis=1)

    if verbose and 'padj' in results_df.columns:
        print('\nThese are all FDR hits:')
        print((results_df['padj'] < fdr).astype(int).groupby('contrast').sum())

    return results_df, full_results_df


def run_linear_models(celltype, model, Y_name, new_coefs, new_coef_dtype, visits, blood, extras, annot_fn, annot_df,
                      visit_interaction, visit_interaction_and_correct_evening, visit_interaction_and_correct_season,
                      LMM, scale, donor_fixed_effects, variance_partition, correct_blood_online,
                      remove_samples, fold_changes, steady_state_season_hack, permutation,
                      fdr, return_corrected_X, save_outputs, save_models, fn_suffix, no_covars=False, verbose=True):

    if isinstance(new_coefs, str):
        new_coefs = [new_coefs]
    if new_coefs is None:
        new_coefs = []

    if variance_partition:
        print('RUNNING VARIANCE PARTITIONING')
        assert LMM and not donor_fixed_effects and not no_covars
        assert len(new_coefs) == 0
    else:
        assert len(new_coefs) != 0
        assert not (fold_changes and LMM)
        assert not (fold_changes and visit_interaction)
        assert not (fold_changes and visit_interaction_and_correct_evening)
        assert not (fold_changes and visit_interaction_and_correct_season)
        assert not (visit_interaction and visit_interaction_and_correct_evening)
        assert not (visit_interaction and visit_interaction_and_correct_season)
        assert not (visit_interaction_and_correct_evening and visit_interaction_and_correct_season)
        assert not ((fold_changes or visit_interaction or visit_interaction_and_correct_evening or visit_interaction_and_correct_season) and any([c.startswith('SAMPLE:') for c in new_coefs if c != 'SAMPLE:VISIT']))

    if verbose and extras:
        print('extras:', extras)

    lmm_groups = None
    if donor_fixed_effects:
        baseline_cols = ['SAMPLE:DONOR', 'SAMPLE:VISIT_TIME_REAL'] + (list(blood) if blood else [])
    else:
        if variance_partition:
            baseline_cols = list(extras) if extras else []
            print('Baseline cols for variance partitioning:', baseline_cols)
        else:
            baseline_cols = sorted(set(BASELINE_COLS[Y_name] + (list(blood) if blood else []) + (list(extras) if extras else [])))
            if not _is_single_visit(visits) and not variance_partition:
                baseline_cols.append('SAMPLE:VISIT')

        if _is_single_visit(visits):
            pass
        elif LMM:
            baseline_cols.append('SAMPLE:DONOR')
            lmm_groups = ['SAMPLE:DONOR']
        elif fold_changes or all([not (c.startswith('DONOR:') or c.endswith('_responder')) for c in new_coefs]):
            baseline_cols.append('SAMPLE:DONOR')
            for c in [b for b in BASELINE_COLS[Y_name] if b.startswith('DONOR:')]:
                baseline_cols.remove(c)
        else:
            raise ValueError

        for c in new_coefs:
            if c in baseline_cols:
                baseline_cols.remove(c)

    if no_covars:
        print(baseline_cols)
        baseline_cols = [c for c in baseline_cols if c in ['SAMPLE:VISIT', 'SAMPLE:DONOR']]

    if 'DONOR:SEX' in baseline_cols and any([c in ['DONOR:postMenopausal', 'DONOR:postMenopausal_binary', 'DONOR:oralContraceptives'] for c in new_coefs]):
        baseline_cols.remove('DONOR:SEX')
    elif 'DONOR:SEX' not in baseline_cols and any(['IncludingMen' in c for c in new_coefs]):
        baseline_cols.append('DONOR:SEX')

    # use all cols for both models to have the same number of vars
    if correct_blood_online and blood:
        correct_columns = ['SAMPLE:VISIT', 'SAMPLE:DONOR'] + blood
        correct_transform = lambda x: _correct_Y(
            X_df=annot_df[correct_columns], Y_df=x,
            design='1 + {}'.format(' + '.join([utils._encode_coef(c, not is_numeric_dtype(annot_df[c].dtype)) for c in correct_columns])),  #  if c != 'SAMPLE:DONOR'
            lmm_groups=None,  # ['SAMPLE:DONOR']
            do_not_correct=[utils._encode_coef(c, not is_numeric_dtype(annot_df[c].dtype), add_categorical_suffix=True) for c in ['SAMPLE:VISIT', 'SAMPLE:DONOR']])
        _Y_transform = correct_transform if Y_transform[Y_name] is None else lambda x: Y_transform[Y_name](correct_transform(x))
    else:
        _Y_transform = Y_transform[Y_name]

    _report_intercept = (fold_changes and set(new_coefs) == {'SAMPLE:VISIT'})
    # use all cols for both models to have the same number of vars
    shared_args = dict(lmm_groups=lmm_groups, scale=scale, fold_changes=fold_changes, columns=baseline_cols + \
                                                                                 new_coefs + \
                                                                                 (['DONOR:IC_EVENING'] if visit_interaction_and_correct_evening else []) + \
                                                                                 (['DONOR:IC_DATE_2PI_SIN', 'DONOR:IC_DATE_2PI_COS'] if visit_interaction_and_correct_season else []),
                       celltype=celltype, annot_fn=annot_fn, Y_name=Y_name, Y_regex=Y_regex[Y_name],
                       Y_transform=_Y_transform, visits=visits, padj_method='fdr_bh',
                       remove_samples=remove_samples,
                       permutation=permutation,
                       do_not_report_intercept=not _report_intercept,
                       reml=not variance_partition, variance_partition=variance_partition,
                       show_raw_distribution=False, save_outputs=save_outputs, save_models=save_models, fn_suffix=fn_suffix, verbose=verbose)

    base_design = ' + '.join(
        [str(1)] + [utils._encode_coef(c, not is_numeric_dtype(annot_df[c].dtype)) for c in baseline_cols if (not lmm_groups or c not in lmm_groups) and
                    (not fold_changes or c not in [fold_changes[0], 'SAMPLE:DONOR'])]
    )

    if visit_interaction or visit_interaction_and_correct_evening or visit_interaction_and_correct_season:
        if visit_interaction_and_correct_evening:
            base_design += ' + DONOR_IC_EVENING + C(SAMPLE_VISIT):DONOR_IC_EVENING'
        if visit_interaction_and_correct_season:
            base_design += ' + DONOR_IC_DATE_2PI_SIN + DONOR_IC_DATE_2PI_COS + C(SAMPLE_VISIT):DONOR_IC_DATE_2PI_SIN + C(SAMPLE_VISIT):DONOR_IC_DATE_2PI_COS'
        if not steady_state_season_hack:
            base_design += ' + ' + ' + '.join([utils._encode_coef(c, not is_numeric_dtype(annot_df[c].dtype)) for c in new_coefs])
        encoded_new_coefs = ['C(SAMPLE_VISIT):{}'.format(utils._encode_coef(c, not is_numeric_dtype(annot_df[c].dtype))) for c in new_coefs]
        contrasts = [['SAMPLE_VISIT\[T\.{visit}\]:{coef}'.format(visit=v, coef=utils._encode_coef(c, not is_numeric_dtype(annot_df[c].dtype), mark_categorical=False, add_categorical_suffix=True)) for c in new_coefs] for v in ['V2', 'V3']]
        do_not_correct = ['C\(SAMPLE_VISIT\)\[T\.{visit}\]:{coef}'.format(visit=v, coef=utils._encode_coef(c, not is_numeric_dtype(annot_df[c].dtype), add_categorical_suffix=True)) for c in new_coefs for v in ['V2', 'V3']]

        # contrasts += [[utils._encode_coef(c, False) for c in ['SAMPLE:VISIT_DATE_2PI_SIN', 'SAMPLE:VISIT_DATE_2PI_COS']], [utils._encode_coef(c, not is_numeric_dtype(annot_df[c].dtype)) for c in new_coefs]]
        # do_not_correct = [utils._encode_coef(c, False) for c in ['SAMPLE:VISIT_DATE_2PI_SIN', 'SAMPLE:VISIT_DATE_2PI_COS']]

        # SOME HACKS TO STUDY SEASONALITY #
        if steady_state_season_hack and set(new_coefs) == set(['DONOR:IC_DATE_2PI_SIN', 'DONOR:IC_DATE_2PI_COS']):
            print('Hacking season analysis')
            _season_coefs = [utils._encode_coef(c, not is_numeric_dtype(annot_df[c].dtype)) for c in ['SAMPLE:VISIT_DATE_2PI_SIN', 'SAMPLE:VISIT_DATE_2PI_COS']]
            contrasts.append(_season_coefs)
            do_not_correct = _season_coefs
            print('contrasts:', contrasts)
            print('do_not_correct:', do_not_correct)
        ##
    elif _report_intercept:
        encoded_new_coefs = []
        contrasts = 'Intercept'
        do_not_correct = list(contrasts)
    else:
        encoded_new_coefs = [utils._encode_coef(c, not is_numeric_dtype(annot_df[c].dtype)) for c in new_coefs]
        contrasts = [utils._encode_coef(c, not is_numeric_dtype(annot_df[c].dtype), mark_categorical=False, add_categorical_suffix=True) for c in new_coefs]
        do_not_correct = [utils._encode_coef(c, not is_numeric_dtype(annot_df[c].dtype), add_categorical_suffix=True) for c in new_coefs]
    full_design = base_design + (' + ' if len(encoded_new_coefs) != 0 else '') + ' + '.join(encoded_new_coefs)

    if verbose:
        for design, design_name in [(base_design, 'base_design'), (full_design, 'full_design')]:
            for _blood in BLOOD_NAMES:
                if len(BLOOD_NAMES[_blood]) != 0:
                    design = design.replace(' + '.join([utils.replace_chars_for_LM(b) for b in BLOOD_NAMES[_blood]]), _blood)
            print('{}:\n~ {}{}'.format(
                design_name, design, ' | groups: {}'.format(lmm_groups) if lmm_groups
                else ' | fold-change: {}'.format(
                    '{}/{} {}'.format(fold_changes[2], fold_changes[1], fold_changes[0])) if fold_changes
                else ''
            ))

    _LR_results = misc.LR_analysis(model=model, design=full_design, contrasts=contrasts,
                                   return_corrected_X=return_corrected_X, do_not_correct=do_not_correct,
                                   var_groups=['PBMC_PERC', 'WB_PER_ML', 'SAMPLE_VISIT', 'SAMPLE_DONOR'] if variance_partition else None,
                                   **shared_args)
    models, results_df, contrasts_df, corrected_df = (_LR_results[0], _LR_results[1], _LR_results[2], _LR_results[3]) if return_corrected_X else (_LR_results[0], _LR_results[1], _LR_results[2], None)

    full_results_df = results_df.copy()

    if len(new_coefs) > 1 or len(new_coefs) == 0:
        results_df = contrasts_df
    else:
        selected = set([])
        for _interaction in ['V2', 'V3'] if visit_interaction or visit_interaction_and_correct_evening or visit_interaction_and_correct_season else [None]:
            select = '^{interaction}{coef}{category}$'.format(
                interaction='SAMPLE_VISIT\[T\.{}\]:'.format(_interaction) if _interaction else '',
                # coef=utils.replace_chars_for_OLS(new_coef) if not _report_intercept else 'Intercept',
                coef=('(?:' + '|'.join(['(?:\\b' + utils.replace_chars_for_LM(c) + '\\b)' for c in new_coefs]) + ')') if not _report_intercept else 'Intercept',
                category='' if new_coef_dtype is None or is_numeric_dtype(new_coef_dtype) or _report_intercept else '\[.*\]'
            )
            mask = results_df.index.get_level_values('contrast').str.contains(select)
            if mask.sum() == 0:
                raise ValueError('Coefficients not found.\nAll contrasts: {}'.format(set(results_df.index.get_level_values('contrast'))))
            else:
                selected = selected.union(results_df.loc[mask].index.get_level_values('contrast'))
        results_df = results_df.loc[results_df.index.get_level_values('contrast').isin(selected)]
    return results_df, full_results_df


if __name__ == '__main__':
    main()
