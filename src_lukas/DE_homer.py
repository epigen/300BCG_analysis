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
    parser.add_argument('--directions', nargs='+', choices=[1, -1], type=int)
    parser.add_argument('--results_dir')
    parser.add_argument('--peak_fn')
    parser.add_argument('--rank_metric', choices=['p.value', 'Coef', 'Coef_times_p.value', 'Coef_times_neg_log10_p.value'])
    parser.add_argument('--effect_size_filter', type=float)
    parser.add_argument('--top_n', type=float, nargs='+')
    parser.add_argument('--metric')
    parser.add_argument('--size', type=int)

    parser.add_argument('--seed', type=int)
    parser.add_argument('--random', action='store_true')
    parser.add_argument('--universe')
    parser.add_argument('--n_jobs', type=int, default=8)

    args = parser.parse_args()

    assert args.contrasts is not None and len(args.contrasts) != 0
    assert args.directions is not None and len(args.directions) != 0
    if args.random and args.seed is None:
        args.seed = int(os.environ.get('SLURM_ARRAY_TASK_ID'))
    assert args.universe is None or os.path.exists(args.universe), args.universe
    assert not args.random or args.seed is not None
    assert not (args.config and (
            args.model or args.celltype
            or args.metric or args.results_dir or args.peak_fn or args.rank_metric or args.effect_size_filter
            or args.top_n))

    if args.config:
         with open(args.config, 'r') as f:
            config = yaml.load(f, Loader=yaml.loader.FullLoader)
            args.model = config['model']
            args.celltype = config['celltype']
            args.results_dir = config['results_dir']
            args.peak_fn = config['peak_fn']
            args.rank_metric = config['rank_metric']
            args.effect_size_filter = config['effect_size_filter']
            args.top_n = config['top_n']
            if not is_iterable(args.top_n):
                args.top_n = [args.top_n]

    if isinstance(args.rank_metric, str):
        args.rank_metric = [args.rank_metric]

    homer_analysis(celltype=args.celltype, model=args.model, contrasts=args.contrasts, directions=args.directions,
                   effect_size_filter=args.effect_size_filter, rank_metrics=args.rank_metric,
                   top_ns=args.top_n, size=args.size, peak_annot_fn=args.peak_fn,
                   universe_fn=args.universe, random=args.random, seed=args.seed,
                   results_dir=args.results_dir, n_jobs=args.n_jobs)


if __name__ == '__main__':
    main()
