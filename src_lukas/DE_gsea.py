import matplotlib as mpl
mpl.use('Agg')
import argparse
import sys
import yaml
import copy
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
    parser.add_argument('--fdr', type=float)
    parser.add_argument('--show_n', type=int)
    parser.add_argument('--metric')
    parser.add_argument('--gene_sets', nargs='+')
    parser.add_argument('--regions', nargs='+', choices=[TSS_PROXIMAL, GENE_AND_DISTAL_10kb, DISTAL_1Mb])
    parser.add_argument('--db_dir')
    parser.add_argument('--nperm', type=int)
    parser.add_argument('--weight', type=float)
    parser.add_argument('--min_size', type=int)
    parser.add_argument('--max_size', type=int)
    parser.add_argument('--nplots', type=float)
    parser.add_argument('--create_svgs', action='store_true')
    parser.add_argument('--random_seed', type=int)
    parser.add_argument('--use_gseapy', action='store_true')
    parser.add_argument('--reduce_to_genes', action='store_true')

    parser.add_argument('--figsize', type=float, nargs='+')
    parser.add_argument('--n_jobs', type=int, default=8)
    parser.add_argument('--memory', default='16G')
    parser.add_argument('--collect_jobs', type=int, default=4)
    parser.add_argument('--collect_memory', default='4G')
    parser.add_argument('--partition', default='shortq')
    parser.add_argument('--walltime', default='12:00:00')
    parser.add_argument('--palette')
    parser.add_argument('--barwidth', type=float, default=0.75)
    parser.add_argument('--plot_context', default='talk', choices=PLOT_CONTEXTS.keys())
    parser.add_argument('--font_scale', type=float, default=1)
    parser.add_argument('--fig_style', default='ticks')
    parser.add_argument('--just_plots', action='store_true')
    parser.add_argument('--run_locally', action='store_true')
    parser.add_argument('--exclude_inodes', action='store_true')
    parser.add_argument('--collect_and_plot', action='store_true')

    args = parser.parse_args()
    assert args.memory[-1].upper() == 'G'
    assert args.collect_memory[-1].upper() == 'G'

    assert not args.figsize or len(args.figsize) == 2
    assert not (args.config and (
            args.model or args.celltype or args.contrasts or args.regions or args.gene_sets or args.show_n
            or args.metric or args.results_dir or args.peak_fn or args.rank_metric or args.effect_size_filter
            or args.db_dir or args.nperm or args.weight or args.min_size or args.max_size or args.fdr or args.nplots
            or args.create_svgs or args.random_seed or args.use_gseapy))

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
            args.show_n = config['gsea_show_n']
            args.metric = config['gsea_metric']
            args.reduce_to_genes = config['gsea_reduce_to_genes']

            args.gene_sets = config['gsea_gene_sets']
            if not args.reduce_to_genes:
                args.gene_sets = ['{}_regions'.format(gs) for gs in
                                  ([args.gene_sets] if isinstance(args.gene_sets, str) else args.gene_sets)]
                print(args.gene_sets)

            args.regions = config['gsea_regions']
            args.fdr = config['gsea_fdr']
            args.db_dir = config['gsea_db_dir']
            args.nperm = config['gsea_nperm']
            args.weight = config['gsea_weight']
            args.min_size = config['gsea_min_size']
            args.max_size = config['gsea_max_size']
            args.nplots = config['gsea_nplots']
            args.create_svgs = config['gsea_create_svgs']
            args.random_seed = config['gsea_random_seed']
            args.use_gseapy = config['gsea_use_gseapy']


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

    if not args.use_gseapy:
        assert args.nplots is None or args.nplots >= 1
    if args.nplots and args.nplots >= 1:
        args.nplots = int(args.nplots)

    if not args.collect_and_plot and not args.just_plots:

        if not args.use_gseapy:
            _args = copy.copy(args)
            _args.config = None
            _args.collect_and_plot = True
            _args.regions = '{regions}'
            unparser = ArgumentUnparser()
            collect_and_plot_cmd = '{} {}'.format(sys.argv[0], unparser.unparse(**vars(_args)))
        else:
            collect_and_plot_cmd = None

        gsea_jobs(celltype=args.celltype, model=args.model, contrasts=args.contrasts,
                  effect_size_filter=args.effect_size_filter, rank_metrics=args.rank_metric, region_filters=args.regions,
                  peak_annot_fn=args.peak_fn, gene_sets=args.gene_sets, db_dir=args.db_dir,
                  nperm=args.nperm, weight=args.weight, min_size=args.min_size, max_size=args.max_size,
                  nplots=args.nplots, create_svgs=args.create_svgs, random_seed=args.random_seed,
                  results_dir=args.results_dir, reduce_to_genes=args.reduce_to_genes,
                  use_gseapy=args.use_gseapy, collect_and_plot_cmd=collect_and_plot_cmd,
                  n_jobs=args.n_jobs, memory=args.memory,
                  collect_jobs=args.collect_jobs, collect_memory=args.collect_memory,
                  partition=args.partition, walltime=args.walltime, exclude_inodes=args.exclude_inodes,
                  run_locally=args.run_locally, show_n=args.show_n, metric=args.metric,
                  palette=args.palette, barwidth=args.barwidth, figsize=args.figsize,
                  gene_sets_fixes=GENE_SETS_FIXES, fdr=args.fdr)

    if args.collect_and_plot and not args.just_plots and not args.use_gseapy:
        for rank_metric in args.rank_metric:
            for region_filter in args.regions:
                for coef in args.contrasts:
                    print('Collect GSEA {} {}'.format(region_filter, coef))
                    gsea_collect(celltype=args.celltype, model=args.model, coef=coef,
                                 effect_size_filter=args.effect_size_filter, rank_metric=rank_metric,
                                 region_filter=region_filter, gene_sets=args.gene_sets, db_dir=args.db_dir,
                                 results_dir=args.results_dir)

    if args.collect_and_plot or args.just_plots or args.use_gseapy:
        fig_axs = plot_gsea(celltype=args.celltype, model=args.model, contrasts=args.contrasts,
                            effect_size_filter=args.effect_size_filter, rank_metrics=args.rank_metric,
                     region_filters=args.regions, gene_sets=None, show_n=args.show_n, metric=args.metric,
                     results_dir=args.results_dir, palette=args.palette, barwidth=args.barwidth, figsize=args.figsize,
                     gene_sets_fixes=GENE_SETS_FIXES, fdr=args.fdr, min_pval=(1 / args.nperm), use_gseapy=args.use_gseapy)
        if fig_axs:
            for p in fig_axs[0]:
                plt.close(fig_axs[0][p])


if __name__ == '__main__':
    main()
