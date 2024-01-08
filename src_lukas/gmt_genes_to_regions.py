import argparse
from misc import get_peak_annot, map_gene_sets_to_regions
from bcg_utils import gene_set_library, read_gene_sets


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--gmt_in')
    parser.add_argument('--gmt_out')
    parser.add_argument('--peak_annot')
    parser.add_argument('--min_size', type=int)
    parser.add_argument('--max_size', type=int)
    parser.add_argument('--save_filtered_gmt')
    parser.add_argument('--NA_symbol', default='--')
    args = parser.parse_args()
    assert args.gmt_in != args.gmt_out

    regions_library = map_gene_sets_to_regions(args.gmt_in,
                                               get_peak_annot(fn=args.peak_annot)['gene_name'].str.upper(),
                                               make_upper=True)

    if args.min_size or args.max_size:
        sizes = {}
        with open(args.save_filtered_gmt if args.save_filtered_gmt else '/dev/null', 'w') as f:
            for term, genes in read_gene_sets(args.gmt_in):
                sizes[term] = len(genes)
                if args.save_filtered_gmt \
                        and (args.min_size is None or sizes[term] >= args.min_size) \
                        and (args.max_size is None or sizes[term] <= args.max_size):
                    print('\t'.join([term, args.NA_symbol] + genes.tolist()), file=f)

    with open(args.gmt_out, 'w') as f:
        for term in regions_library:
            if (args.min_size is None or sizes[term] >= args.min_size) \
                    and (args.max_size is None or sizes[term] <= args.max_size):
                print('\t'.join([term, args.NA_symbol] + regions_library[term]), file=f)

    # Test if it is still in a correct format:
    gene_set_library(fn=args.gmt_out, description_NA=args.NA_symbol)
    if args.save_filtered_gmt:
        gene_set_library(fn=args.save_filtered_gmt, description_NA=args.NA_symbol)


if __name__ == '__main__':
    main()
