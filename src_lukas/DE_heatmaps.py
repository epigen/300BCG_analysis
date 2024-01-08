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
    parser.add_argument('--do_not_correct', nargs='+')
    parser.add_argument('--annotations', nargs='+')
    parser.add_argument('--palettes', nargs='+', choices=PALETTE_CODES.keys())
    parser.add_argument('--outfile')
    parser.add_argument('--annot_fn')
    parser.add_argument('--rank_metric', choices=['p.value', 'Coef'])
    parser.add_argument('--top_n', type=float, nargs='+')

    parser.add_argument('--figsize', type=float, nargs='+', default=[8, 8])
    parser.add_argument('--legend_bbox', type=float, nargs='+')
    parser.add_argument('--figsize_counts', type=float, nargs='+', default=[5, 4.5])
    parser.add_argument('--dpi', type=int, default=300)
    parser.add_argument('--not_rasterized', action='store_true')
    parser.add_argument('--plot_context', choices=PLOT_CONTEXTS.keys(), default='talk')
    parser.add_argument('--font_scale', type=float, default=1)
    parser.add_argument('--fig_style', default='ticks')
    args = parser.parse_args()

    assert not args.figsize or len(args.figsize) == 2
    assert not args.figsize_counts or len(args.figsize_counts) == 2
    assert not (args.config and (
            args.model or args.celltype or args.contrasts \
            or args.results_dir or args.outfile or args.do_not_correct \
            or args.annot_fn or args.top_n or args.annotations or args.palettes))

    if args.config:
         with open(args.config, 'r') as f:
            config = yaml.load(f, Loader=yaml.loader.FullLoader)
            args.model = config['model']
            args.celltype = config['celltype']
            args.contrasts = config['contrasts']
            args.results_dir = config['results_dir']
            args.do_not_correct = config['do_not_correct']
            args.rank_metric = config['rank_metric']
            args.top_n = config['top_n']
            if not is_iterable(args.top_n):
                args.top_n = [args.top_n]
            args.annot_fn = config['annot_fn']
            args.annotations = config['heatmap_annotations']
            args.palettes = config['heatmap_palettes']

    if isinstance(args.contrasts, str):
        args.contrasts = [args.contrasts]
    if is_iterable(args.contrasts[0]):
        assert len(args.contrasts) == 1
        args.contrasts = args.contrasts[0]
    if isinstance(args.top_n, int) or isinstance(args.top_n, float):
        args.top_n = [args.top_n]

    if args.outfile is None:
        args.outfile = de_fn(celltype=args.celltype, model=args.model,
                             data='{pseudo_sort}_{plot}',
                             extras='{coef}.{rank_metric}.{top_n}.{direction}',
                             ext='pdf', gzipped=False, subdir='heatmaps', results_dir=args.results_dir)

    setup_plotting(style=args.fig_style, context=args.plot_context, font_scale=args.font_scale)

    atac_df, _ = get_corrected_counts(celltype=args.celltype, model=args.model, do_not_correct=args.do_not_correct,
                                      remove_singletons=False, plot=True, figsize=args.figsize_counts,
                                      fig_fn=de_fn(celltype=args.celltype, model=args.model,
                                                   data='normalized_corrected_counts',
                                                   ext='pdf', gzipped=False, results_dir=args.results_dir))
    plt.close()

    annot_df = pd.read_csv(args.annot_fn, index_col=0).loc[atac_df.index]
    de_df = read_de(celltype=args.celltype, model=args.model, contrasts=args.contrasts)

    for coef in args.contrasts:
        for top_n in args.top_n:
            for pseudo_sort in [True, False]:
                if not pseudo_sort or coef in ['V3', 'V2']:
                    fig_fn = args.outfile.format(
                        coef=coef, rank_metric=args.rank_metric,
                        top_n='{}{}'.format('top_' if top_n >= 1 else 'FDR_', top_n),
                        direction='{direction}',
                        pseudo_sort='pseudosort' if pseudo_sort else 'cluster',
                        plot='{plot}'
                    )
                    try:
                        cg = de_heatmaps(atac_df, annot_df, de_df, coef, args.rank_metric, top_n,
                                         args.annotations, palettes=args.palettes,
                                         pseudo_sort=pseudo_sort, pseudosort_col='SAMPLE:VISIT', pseudosort_base_value='V1',
                                         pseudosort_hist_height=0.75,
                                         pseudosort_hist_fn=fig_fn.format(direction='{direction}', plot='pseudosort_hist'),
                                         rasterized=not args.not_rasterized, legend_bbox=args.legend_bbox if args.legend_bbox else 'auto',
                                         figsize=args.figsize)
                        if cg:
                            if not pseudo_sort:
                                cg = {'Up/Down': cg}
                            for direction in cg:
                                cg[direction].ax_col_colors.tick_params(labelleft=True, labelright=False, labelrotation=0)
                                if args.annotations and len(args.annotations) == 1:
                                    # cg[direction].ax_col_colors.set_yticklabels([])
                                    pass
                                savefig(fig_fn.format(direction=direction.lower().replace('/', '_'),
                                                      plot='{}_heatmap'.format('pseudosort' if pseudo_sort else 'cluster')),
                                        fig=cg[direction].fig, dpi=args.dpi)
                                plt.close(cg[direction].fig)
                    except ValueError as e:
                        plt.close('all')
                        print('ERROR ({} {} {}): {}'.format(coef, top_n, 'pseudosort' if pseudo_sort else 'cluster', e))


if __name__ == '__main__':
    main()
