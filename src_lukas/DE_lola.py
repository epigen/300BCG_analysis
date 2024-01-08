import matplotlib as mpl
mpl.use('Agg')
import argparse
import yaml
from misc import *
from bcg_utils import *


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--config')
    parser.add_argument('--model')
    parser.add_argument('--celltype')
    parser.add_argument('--contrasts', nargs='+')
    parser.add_argument('--results_dir')
    parser.add_argument('--peak_fn')
    parser.add_argument('--rank_metric', choices=['p.value', 'Coef', 'Coef_times_p.value', 'Coef_times_neg_log10_p.value'])
    parser.add_argument('--effect_size_filter', type=float)
    parser.add_argument('--top_n', type=float, nargs='+')
    parser.add_argument('--fdr', type=float)
    parser.add_argument('--show_n', type=int)
    parser.add_argument('--metric')
    parser.add_argument('--db')
    parser.add_argument('--collections')
    parser.add_argument('--remove_replicates', action='store_true')
    parser.add_argument('--regions', nargs='+', choices=[TSS_PROXIMAL, GENE_AND_DISTAL_10kb, DISTAL_1Mb, PROMOTERS, ALL_GENES])
    parser.add_argument('--db_dir')

    parser.add_argument('--seed', type=int)
    parser.add_argument('--random', action='store_true')
    parser.add_argument('--universe')

    parser.add_argument('--figsize', type=float, nargs='+')
    parser.add_argument('--figsize_distributions', type=float, nargs='+', default=[5, 5])
    parser.add_argument('--n_jobs', type=int, default=8)
    parser.add_argument('--palette')
    parser.add_argument('--barwidth', type=float, default=0.75)
    parser.add_argument('--plot_context', default='talk', choices=PLOT_CONTEXTS.keys())
    parser.add_argument('--font_scale', type=float, default=1)
    parser.add_argument('--fig_style', default='ticks')
    parser.add_argument('--just_plots', action='store_true')
    parser.add_argument('--skip_lola', action='store_true')
    args = parser.parse_args()

    if args.random and args.seed is None:
        args.seed = int(os.environ.get('SLURM_ARRAY_TASK_ID'))
    assert args.universe is None or os.path.exists(args.universe), args.universe
    assert not args.random or args.seed is not None
    assert not args.figsize or len(args.figsize) == 2
    assert not args.figsize_distributions or len(args.figsize_distributions) == 2
    assert not (args.config and (
            args.model or args.celltype or args.contrasts or args.db or args.collections or args.show_n
            or args.metric or args.results_dir or args.peak_fn or args.rank_metric or args.effect_size_filter
            or args.db_dir or args.remove_replicates or args.fdr))

    if args.config:
         with open(args.config, 'r') as f:
            config = yaml.load(f, Loader=yaml.loader.FullLoader)
            args.model = config['model']
            args.celltype = config['celltype']
            args.contrasts = config['enr_contrasts'] if 'enr_contrasts' in config else config['contrasts']
            args.results_dir = config['results_dir']
            args.peak_fn = config['peak_fn']
            args.rank_metric = config['rank_metric']
            args.effect_size_filter = config['effect_size_filter']
            if args.top_n is None:
                args.top_n = config['top_n']
                if not is_iterable(args.top_n):
                    args.top_n = [args.top_n]
            args.db = config['lola_db']
            args.show_n = config['lola_show_n']
            args.metric = config['lola_metric']
            args.collections = config['lola_collections']
            args.db_dir = config['lola_db_dir']
            args.remove_replicates = config['lola_remove_replicates']
            args.fdr = config['lola_fdr']

    if isinstance(args.contrasts, str):
        args.contrasts = [args.contrasts]
    if isinstance(args.rank_metric, str):
        args.rank_metric = [args.rank_metric]
    if is_iterable(args.contrasts[0]):
        assert len(args.contrasts) == 1
        args.contrasts = args.contrasts[0]
    if args.figsize is None:
        args.figsize = [10, 0.4 * args.show_n]

    setup_plotting(style=args.fig_style, context=args.plot_context, font_scale=args.font_scale)
    
    if not args.just_plots:
        lola_analysis(celltype=args.celltype, model=args.model, contrasts=args.contrasts, region_filters=args.regions,
                      effect_size_filter=args.effect_size_filter, rank_metrics=args.rank_metric, top_ns=[int(top_n) if top_n >= 1 else top_n for top_n in args.top_n],
                      peak_annot_fn=args.peak_fn, universe_fn=args.universe, databases=args.db, db_dir=args.db_dir,
                      run_lola=(not args.skip_lola), results_dir=args.results_dir,
                      palette=args.palette if args.palette else 'deep', figsize=args.figsize_distributions,
                      plot=not args.random, random=args.random, seed=args.seed, n_jobs=args.n_jobs)

#    if False and not args.skip_lola and not args.random:
#        plot_lola(celltype=args.celltype, model=args.model, contrasts=args.contrasts,
#                  effect_size_filter=args.effect_size_filter, rank_metrics=args.rank_metric, top_ns=[int(top_n) if top_n >= 1 else top_n for top_n in args.top_n],
#                  databases=None, collections=args.collections, show_n=args.show_n,
#                  metric=args.metric, remove_replicates=args.remove_replicates, results_dir=args.results_dir,
#                  palette=args.palette, barwidth=args.barwidth, figsize=args.figsize,
#                  celltype_fixes=LOLA_CELLTYPES_FIXES, fdr=args.fdr)


if __name__ == '__main__':
    main()
