import argparse
from bcg_utils import read_gene_sets, gene_set_library


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--fn')
    parser.add_argument('--out')
    parser.add_argument('--NA_symbol', default='--')
    parser.add_argument('--make_safe_names', action='store_true', help='This will replace "/" with "-"')
    args = parser.parse_args()
    assert args.fn != args.out
    assert len(args.NA_symbol) != 0

    _enrichr_to_gmt(fn=args.fn, out=args.out, NA_symbol=args.NA_symbol, make_safe_names=args.make_safe_names)
    # Test if it is still in a correct format:
    gene_set_library(fn=args.out, description_NA=args.NA_symbol)


def _enrichr_to_gmt(fn, out, NA_symbol, make_safe_names):
    with open(out, 'w') as f:
        for term, genes in read_gene_sets(fn=fn, description_NA=None, as_arrays=False, assert_upper=True):
            if make_safe_names and '/' in term:
                print('Fixing: "{}"'.format(term))
                term = term.replace('/', '-')
                print('Fixed: "{}"'.format(term))
            print('\t'.join([term, NA_symbol] + genes), file=f)


if __name__ == '__main__':
    main()
