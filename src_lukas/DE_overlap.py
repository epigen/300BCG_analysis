import matplotlib as mpl
mpl.use('Agg')
import argparse
import yaml
from misc import *
# from bcg_utils import *


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--config')
    parser.add_argument('--model')
    parser.add_argument('--compare_celltype')
    parser.add_argument('--compare_model')
    parser.add_argument('--compare_legacy', action='store_true')
    parser.add_argument('--celltype')
    parser.add_argument('--contrasts', nargs='+')
    parser.add_argument('--compare_contrasts', nargs='+')
    parser.add_argument('--results_dir')
    parser.add_argument('--outfile')

    parser.add_argument('--figsize', type=float, nargs='+', default=[5, 4.5])
    parser.add_argument('--not_rasterized', action='store_true')
    parser.add_argument('--plot_context', choices=PLOT_CONTEXTS.keys(), default='talk')
    parser.add_argument('--font_scale', type=float, default=1)
    parser.add_argument('--fig_style', default='ticks')
    args = parser.parse_args()

    assert not args.figsize or len(args.figsize) == 2
    assert args.compare_contrasts is None or len(args.compare_contrasts) == len(args.contrasts)
    assert not (args.config and (
            args.model or args.compare_model or args.compare_celltype or args.compare_legacy or args.celltype \
            or args.contrasts or args.compare_contrasts or args.results_dir or args.outfile))

    if args.config:
         with open(args.config, 'r') as f:
            config = yaml.load(f, Loader=yaml.loader.FullLoader)
            args.model = config['model']
            args.compare_model = config['compare_model']
            args.compare_celltype = config['compare_celltype']
            args.compare_legacy = config['compare_legacy']
            args.celltype = config['celltype']
            args.contrasts = config['contrasts']
            args.compare_contrasts = config['compare_contrasts']
            args.results_dir = config['results_dir']

    if isinstance(args.contrasts, str):
        args.contrasts = [args.contrasts]
    if is_iterable(args.contrasts[0]):
        assert len(args.contrasts) == 1
        args.contrasts = args.contrasts[0]

    if args.outfile is None:
        args.outfile = de_fn(celltype=args.celltype, model=args.model, data='{data}',
                             ext='pdf', gzipped=False, results_dir='results')

    if args.compare_contrasts is None:
        args.compare_contrasts = args.contrasts

    setup_plotting(style=args.fig_style, context=args.plot_context, font_scale=args.font_scale)
    df = read_de(celltype=args.celltype, model=args.model, contrasts=args.contrasts)
    baseline_df = read_de(celltype=args.compare_celltype if args.compare_celltype else args.celltype,
                          model=args.compare_model, contrasts=args.compare_contrasts, legacy=args.compare_legacy)

    _ = plot_de_agreement(df, baseline_df, args.contrasts, args.compare_contrasts, figsize=args.figsize,
                          rasterized=not args.not_rasterized)
    savefig(args.outfile.format(data='agreement'))
    plt.close()

    overlap = {}
    for coef, baseline_coef in zip(args.contrasts, args.compare_contrasts):
        overlap[coef] = de_overlap(df, baseline_df, coef, baseline_coef=baseline_coef)
    _ = plot_de_overlap(overlap, figsize=args.figsize)
    savefig(args.outfile.format(data='overlap'))
    plt.close()


if __name__ == '__main__':
    main()
