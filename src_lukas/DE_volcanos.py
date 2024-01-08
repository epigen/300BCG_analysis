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
    parser.add_argument('--hue', default='A')
    parser.add_argument('--fdr', type=float)
    parser.add_argument('--outfile')

    parser.add_argument('--figsize', type=float, nargs='+', default=[5, 4.5])
    parser.add_argument('--dpi', type=int, default=300)
    parser.add_argument('--sharex', action='store_true')
    parser.add_argument('--sharey', action='store_true')
    parser.add_argument('--not_rasterized', action='store_true')
    parser.add_argument('--plot_context', default='talk', choices=PLOT_CONTEXTS.keys())
    parser.add_argument('--font_scale', type=float, default=1)
    parser.add_argument('--fig_style', default='ticks')

    args = parser.parse_args()
    assert not args.figsize or len(args.figsize) == 2
    assert not (args.config and (
            args.model or args.celltype or args.contrasts \
            or args.results_dir or args.outfile or args.fdr))

    if args.config:
         with open(args.config, 'r') as f:
            config = yaml.load(f, Loader=yaml.loader.FullLoader)
            args.model = config['model']
            args.celltype = config['celltype']
            args.contrasts = config['contrasts']
            args.results_dir = config['results_dir']
            args.fdr = config['volcano_fdr']

    if isinstance(args.contrasts, str):
        args.contrasts = [args.contrasts]

    setup_plotting(style=args.fig_style, context=args.plot_context, font_scale=args.font_scale)

    if 'dream' in args.model.lower():
        assert len(args.contrasts) == 1
        _first = True
        de_df = []
        for contrast in args.contrasts[0]:
            fn = de_fn(args.celltype, args.model, data='results_{}'.format(contrast))
            df = pd.read_csv(fn, index_col=0)
            df.columns = df.columns.str.replace('^logFC$', 'Coef.{}'.format(contrast)).str.replace(
                '^AveExpr$', 'A').str.replace('^t$', 't.{}'.format(contrast)).str.replace(
                '^P.Value$', 'p.value.{}'.format(contrast)).str.replace('^adj.P.Val$', 'padj.{}'.format(contrast))
            df = df.drop(['z.std', 'genes'], axis=1)
            if not _first:
                df = df.drop(['A'], axis=1)
            de_df.append(df)
            _first = False
        de_df = pd.concat(de_df, axis=1)
        de_df.to_csv(de_fn(args.celltype, args.model, data='results_p5'))

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

        df = read_de(celltype=args.celltype, model=args.model, contrasts=contrasts,
                     results_dir=args.results_dir, extras=_contrasts_suffix, add_coef_times_pval=True)
        axs = volcano(df, contrasts=contrasts, hue=args.hue, fdr=args.fdr, figsize=args.figsize, rasterized=not args.not_rasterized,
                      sharex=args.sharex, sharey=args.sharey,
                      legend_kwargs=dict(bbox_to_anchor=(1, 1),
                                         title='Average\nlog$_2$CPM' if args.hue == 'A' else args.hue))
        for col, coef in enumerate(contrasts):
            if coef in ['V2', 'V3']:
                try:
                    effective_n = get_de_effective_N(args.celltype, args.model, visit1='V1', visit2=coef)
                    axs[col].set_title('{} ($n$={})'.format(coef, effective_n))
                except FileNotFoundError:
                    pass
        savefig(args.outfile if args.outfile else de_fn(
            celltype=args.celltype, model=args.model, data='volcano', extras=_contrasts_suffix,
            ext='pdf', gzipped=False, results_dir=args.results_dir
        ), dpi=args.dpi)
        plt.close()


if __name__ == '__main__':
    main()
