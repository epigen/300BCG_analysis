import argparse
import bcg_utils as utils


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--data_fn', required=True)
    parser.add_argument('--celltype_col', required=True)
    parser.add_argument('--reference')
    parser.add_argument('--interaction_model', action='store_true')
    parser.add_argument('--subsample_celltypes', action='store_true')
    parser.add_argument('--only_epi_genes', action='store_true')
    parser.add_argument('--results_dir', required=True)
    args = parser.parse_args()

    utils.sclm2(data_fn=args.data_fn, celltype_col=args.celltype_col,
               results_dir=args.results_dir, reference=args.reference,
               interaction_model=args.interaction_model, subsample_celltypes=args.subsample_celltypes,
                only_epi_genes=args.only_epi_genes)


if __name__ == '__main__':
    main()
