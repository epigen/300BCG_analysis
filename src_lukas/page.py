import os
import pandas as pd
import numpy as np
from scipy.stats import norm, rankdata
from statsmodels.stats.multitest import multipletests


def _read_gene_sets(fn, NA_symbol='--', as_arrays=True, assert_upper=True):
    '''
    Reads in a gene set file in a GMT format.

    :param fn: gene set filename downloaded from Enrichr website
    :param NA_symbol: I assume description is missing or equal to "N/A" symbol for safety
    :param as_arrays: whether the second item of the tuple should be an array or a list
    :param assert_upper: fail if there is a gene with lower-case characters
    :return: iterable of tuples (term, genes)
    '''

    with open(fn) as f:
        for line in f:
            line_split = line.strip().split('\t')
            term, description, genes = line_split[0].strip(), line_split[1].strip(), line_split[2:]
            assert len(term) != 0
            assert len(description) == 0 or description == NA_symbol
            assert all([' ' not in g and ',' not in g and (not assert_upper or g == g.upper()) for g in genes])
            yield term, np.asarray(genes) if as_arrays else genes


def page_test(LFCs, gene_set_mask, min_genes=None, max_genes=None, two_sided=True, ddof=0, axis=1):
    '''
    Performs the PAGE test for one gene set but possible across many comparisons.

    :param LFCs: a dataframe of log2 fold changes; when axis=1, genes are as columns and comparisons as rows
    :param gene_set_mask: a set of genes or a boolean mask for the genes in the LFCs dataframe
    :param min_genes: do not test if less than min_genes would be tested for the given gene set
    :param max_genes: do not test if more than max_genes would be tested for the given gene set
    :param two_sided: two-sided T-test
    :param ddof: degrees of freedom
    :param axis: 1 if genes are in the columns and 0 otherwise
    :return: a tuple of arrays with a value for each comparison (z_score, p_value, mean, std,
    gene_set_mean, gene_set_sizes, gene_set_mask)
    '''

    if axis == 1:
        # genes need to be along the axis 0
        LFCs = LFCs.T

    if gene_set_mask.dtype != bool:
        gene_set_mask = LFCs.index.isin(gene_set_mask)

    # in every comparison a different gene set mask may be applied because of NaNs
    # thus blowing up gene_set_mask 1-D array to a 2-D array
    if (min_genes is None or gene_set_mask.sum() >= min_genes) \
            and (max_genes is None or gene_set_mask.sum() <= max_genes):
        gene_set_mask = np.broadcast_to(gene_set_mask, (LFCs.shape[1], gene_set_mask.shape[0])).T
        gene_set_mask = gene_set_mask & ~np.isnan(LFCs)
        gene_set_sizes = gene_set_mask.sum(axis=0)

        mean = LFCs.mean(axis=0)
        std = np.sqrt(LFCs.var(axis=0, ddof=ddof) / gene_set_sizes)
        gene_set_mean = LFCs[gene_set_mask].mean(axis=0)
        z = (gene_set_mean - mean) / std
        p = norm.sf(np.abs(z))

        return (z, (p * 2) if two_sided else p, mean, std, gene_set_mean,
                gene_set_sizes, gene_set_mask.T if axis == 1 else gene_set_mask)
    else:
        return None, None, None, None, None, gene_set_mask.sum(), None


def page(LFCs, gene_set_libraries, gs_path=None, min_genes=None, max_genes=None, adjust_pvals=True, rank=True):
    '''
    Runs page tests across several gene set libraries and comparisons.

    :param LFCs: a dataframe of log2 fold changes; when axis=1, genes are as columns and comparisons as rows
    :param gene_set_libraries: names of the gene set libraries files
    :param gs_path: directory containing gene set libraries files
    :param min_genes: do not test if less than min_genes would be tested for the given gene set
    :param max_genes: do not test if more than max_genes would be tested for the given gene set
    :param adjust_pvals: perform Benjamini-Hochberg FDR correction
    :param rank: rank by absolute z-score
    :return: a long dataframe of all comparisons, gene set libraries, and terms
    '''

    if isinstance(gene_set_libraries, str):
        gene_set_libraries = [gene_set_libraries]

    page_df = pd.DataFrame()

    for gene_set_library in gene_set_libraries:
        lib_df = pd.DataFrame()
        for term, genes in _read_gene_sets(os.path.join(gs_path, gene_set_library)):
            z_scores, p_values, mean, std, gene_set_mean, gene_set_sizes, gene_set_mask = \
                page_test(LFCs, genes, min_genes=min_genes, max_genes=max_genes)
            if z_scores is not None:
                term_df = pd.DataFrame(index=pd.MultiIndex.from_product(
                    [[gene_set_library], [term], z_scores.index],
                    names=['gene_set_library', 'description', 'comparison']
                ))
                term_df['z_score'] = z_scores.values
                term_df['p_value'] = p_values
                term_df['mean'] = mean.values
                term_df['std'] = std.values
                term_df['gene_set_mean'] = gene_set_mean.values
                term_df['delta_mean'] = (gene_set_mean - mean).values
                term_df['gene_set_size'] = gene_set_sizes.values

                # I was considering defining lead genes as all those that have LFCs larger than the mean
                # but I decided it against it in the end
                # lead_mask = (((LFCs.where(z_scores > 0).T >= np.clip(mean, a_min=0, a_max=None)) | (
                #             LFCs.where(z_scores < 0).T <= np.clip(mean, a_min=None, a_max=0))) & gene_set_mask.T).T

                all_genes, contrib_genes, lead_genes = [], [], []
                for c in z_scores.index:
                    gene_set_LFCs = LFCs.loc[c, gene_set_mask.loc[c]].sort_values(ascending=z_scores.loc[c] < 0)
                    all_genes.append(str(gene_set_LFCs.index.tolist()))
                    contrib_genes.append(str(gene_set_LFCs.loc[(gene_set_LFCs > 0) if z_scores.loc[c] > 0 else (
                                gene_set_LFCs < 0)].index.tolist()))

                    # I was considering defining lead genes as all those that have LFCs larger than the mean
                    # but I decided it against it in the end
                    # lead_LFCs = LFCs.loc[c, lead_mask.loc[c]].sort_values(ascending=z_scores.loc[c] < 0)

                    # lead genes are those with within gene set absolute z-score more than 1
                    # Note: my personal definition
                    lead_LFCs = gene_set_LFCs
                    z_lead_LFCs = (lead_LFCs - lead_LFCs.mean()) / lead_LFCs.std()
                    lead_genes.append(str(lead_LFCs.loc[(z_lead_LFCs >= 1) if z_scores.loc[c] > 0 else (
                                z_lead_LFCs <= -1)].index.tolist()))
                term_df['ledge_genes'] = lead_genes  # leading edge genes
                term_df['contrib_genes'] = contrib_genes  # gene set genes that were up/down if the pathway was up/down
                term_df['genes'] = all_genes  # all the genes from the gene set that were measured
                lib_df = pd.concat([lib_df, term_df], axis=0)

        if adjust_pvals:
            lib_df['adjusted_p_value'] = np.nan
            pval_mask = ~pd.isnull(lib_df['p_value'])
            for comparison in np.unique(lib_df.index.get_level_values('comparison')):
                cmp_mask = lib_df.index.get_level_values('comparison') == comparison
                lib_df.loc[cmp_mask & pval_mask, 'adjusted_p_value'] = \
                    multipletests(lib_df.loc[cmp_mask & pval_mask, 'p_value'].values, method='fdr_bh')[1]

        if rank:
            lib_df['rank'] = np.nan
            for comparison in np.unique(lib_df.index.get_level_values('comparison')):
                cmp_mask = lib_df.index.get_level_values('comparison') == comparison
                lib_df.loc[cmp_mask, 'rank'] = rankdata(-lib_df.loc[cmp_mask, 'z_score'].abs().values, method='ordinal')

        page_df = pd.concat([page_df, lib_df], axis=0)

    return page_df
