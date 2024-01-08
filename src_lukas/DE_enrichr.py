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
    parser.add_argument('--gene_sets', nargs='+')
    parser.add_argument('--regions', nargs='+', choices=[TSS_PROXIMAL, GENE_AND_DISTAL_10kb, DISTAL_1Mb, PROMOTERS, ALL_GENES])
    parser.add_argument('--db_dir')

    parser.add_argument('--figsize', type=float, nargs='+')
    parser.add_argument('--palette')
    parser.add_argument('--barwidth', type=float, default=0.75)
    parser.add_argument('--plot_context', default='talk', choices=PLOT_CONTEXTS.keys())
    parser.add_argument('--font_scale', type=float, default=1)
    parser.add_argument('--fig_style', default='ticks')
    parser.add_argument('--just_plots', action='store_true')
    parser.add_argument('--no_plots', action='store_true')
    parser.add_argument('--F_test', action='store_true')
    args = parser.parse_args()

    assert not args.figsize or len(args.figsize) == 2
    assert not (args.config and (
            args.model or args.celltype or args.contrasts or args.db_dir or args.show_n
            or args.metric or args.results_dir or args.peak_fn or args.rank_metric or args.effect_size_filter # or args.top_n or args.gene_sets or args.regions
            or args.fdr))

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
            args.top_n = config['top_n'] if not args.top_n else args.top_n
            if not is_iterable(args.top_n):
                args.top_n = [args.top_n]
            args.show_n = config['enrichr_show_n']
            args.metric = config['enrichr_metric']
            args.gene_sets = config['enrichr_gene_sets'] if not args.gene_sets else args.gene_sets
            args.db_dir = config['enrichr_db_dir'] if 'enrichr_db_dir' in config else os.path.join(METADATA, 'gene_set_libraries')
            args.regions = config['enrichr_regions'] if not args.regions else args.regions
            args.fdr = config['enrichr_fdr']

    args.top_n = [int(_n) if _n >= 1 else _n for _n in args.top_n]
    if isinstance(args.contrasts, str):
        args.contrasts = [args.contrasts]

    for contrasts in args.contrasts:
        _contrasts_suffix = '_{}'.format('_'.join(contrasts)) if len(args.contrasts) > 1 else None
        if len(contrasts) == 1:
            _fn = de_fn(celltype=args.celltype, model=args.model,
                        results_dir=args.results_dir, extras=_contrasts_suffix)
            df = pd.read_csv(_fn, index_col=0)
            _cols = df.columns.tolist()
            df = df.rename({contrasts[0]: 'Coef.{}'.format(contrasts[0]),
                            '{}.1'.format(contrasts[0]): 't.{}'.format(contrasts[0]),
                            '{}.2'.format(contrasts[0]): 'p.value.{}'.format(contrasts[0]),
                            '{}.3'.format(contrasts[0]): 'Res.{}'.format(contrasts[0])}, axis=1)
            if _cols != df.columns.tolist():
                print('Fixing the columns names in', _fn)
                df.to_csv(_fn)

    if isinstance(args.rank_metric, str):
        args.rank_metric = [args.rank_metric]
    if is_iterable(args.contrasts[0]):
        assert len(args.contrasts) == 1
        args.contrasts = args.contrasts[0]
    if args.figsize is None:
        args.figsize = [14, 0.4 * args.show_n]

    setup_plotting(style=args.fig_style, context=args.plot_context, font_scale=args.font_scale)

    if not args.just_plots:
        # enrichr_online(celltype=args.celltype, model=args.model, contrasts=args.contrasts,
        #               effect_size_filter=args.effect_size_filter, rank_metrics=args.rank_metric, top_ns=args.top_n, region_filters=args.regions,
        #               peak_annot_fn=args.peak_fn, gene_sets=args.gene_sets, results_dir=args.results_dir)

        run_gr_enrichr(celltype=args.celltype, model=args.model, contrasts=args.contrasts,
                   effect_size_filter=args.effect_size_filter, rank_metrics=args.rank_metric, top_ns=args.top_n, region_filters=args.regions,
                   peak_annot_fn=args.peak_fn, gene_sets=args.gene_sets, db_dir=args.db_dir,
                   results_dir=args.results_dir, F_test=args.F_test, strict_background=True, padj_union=False)

    if not args.no_plots:
        fig_axs = plot_enrichr(celltype=args.celltype, model=args.model, contrasts=args.contrasts,
                               effect_size_filter=args.effect_size_filter, rank_metrics=args.rank_metric,
                               top_ns=args.top_n, region_filters=args.regions, gene_sets=None, show_n=args.show_n, metric=args.metric,
                     results_dir=args.results_dir, palette=args.palette, barwidth=args.barwidth, figsize=args.figsize,
                     gene_sets_fixes=GENE_SETS_FIXES, fdr=args.fdr)
        if fig_axs:
            for p in fig_axs[0]:
                plt.close(fig_axs[0][p])


if __name__ == '__main__':
    main()
