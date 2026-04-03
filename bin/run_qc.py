#!/usr/bin/env python3
import argparse
import os
import numpy as np
import scanpy as sc
import scrublet as scr
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt


def run_doublet_detection(adata):
    adata.obs['doublet_score'] = 0.0
    adata.obs['predicted_doublet'] = False

    for sample in adata.obs['sample_id'].unique():
        mask = adata.obs['sample_id'] == sample
        subset = adata[mask].copy()

        scrub = scr.Scrublet(subset.X)
        doublet_scores, predicted_doublets = scrub.scrub_doublets(verbose=False)

        adata.obs.loc[mask, 'doublet_score'] = doublet_scores
        adata.obs.loc[mask, 'predicted_doublet'] = predicted_doublets

        n_doublets = predicted_doublets.sum()
        print(f"  {sample}: {n_doublets}/{mask.sum()} doublets "
              f"({100 * n_doublets / mask.sum():.1f}%)")

    return adata


def run_qc(input_file, output_file, plot_dir, min_genes, min_cells, max_mito_pct):
    os.makedirs(plot_dir, exist_ok=True)

    adata = sc.read_h5ad(input_file)
    print(f"Loaded: {adata.n_obs} cells x {adata.n_vars} genes")

    adata.var['mt'] = adata.var_names.str.startswith('MT-')
    sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], inplace=True)

    fig, axes = plt.subplots(1, 3, figsize=(12, 4))
    sc.pl.violin(adata, 'n_genes_by_counts', ax=axes[0], show=False)
    sc.pl.violin(adata, 'total_counts', ax=axes[1], show=False)
    sc.pl.violin(adata, 'pct_counts_mt', ax=axes[2], show=False)
    plt.tight_layout()
    plt.savefig(os.path.join(plot_dir, 'qc_violin_pre_filter.png'), dpi=150)
    plt.close()

    fig, axes = plt.subplots(1, 2, figsize=(10, 4))
    sc.pl.scatter(adata, x='total_counts', y='n_genes_by_counts', ax=axes[0], show=False)
    sc.pl.scatter(adata, x='total_counts', y='pct_counts_mt', ax=axes[1], show=False)
    plt.tight_layout()
    plt.savefig(os.path.join(plot_dir, 'qc_scatter_pre_filter.png'), dpi=150)
    plt.close()

    n_before = adata.n_obs
    sc.pp.filter_cells(adata, min_genes=min_genes)
    sc.pp.filter_genes(adata, min_cells=min_cells)
    print(f"Filtered: {adata.n_obs} cells remain ({n_before - adata.n_obs} removed, "
          f"min_genes={min_genes}, min_cells={min_cells})")

    n_before = adata.n_obs
    adata = adata[adata.obs.pct_counts_mt < max_mito_pct].copy()
    print(f"Mito filter (<{max_mito_pct}%): {adata.n_obs} cells ({n_before - adata.n_obs} removed)")

    adata = run_doublet_detection(adata)
    n_before = adata.n_obs
    adata = adata[~adata.obs['predicted_doublet']].copy()
    print(f"Doublets removed: {n_before - adata.n_obs}")

    fig, axes = plt.subplots(1, 3, figsize=(12, 4))
    sc.pl.violin(adata, 'n_genes_by_counts', ax=axes[0], show=False)
    sc.pl.violin(adata, 'total_counts', ax=axes[1], show=False)
    sc.pl.violin(adata, 'pct_counts_mt', ax=axes[2], show=False)
    plt.tight_layout()
    plt.savefig(os.path.join(plot_dir, 'qc_violin_post_filter.png'), dpi=150)
    plt.close()

    adata.write_h5ad(output_file)
    print(f"Final: {adata.n_obs} cells x {adata.n_vars} genes -> {output_file}")


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--input', required=True)
    parser.add_argument('--output', required=True)
    parser.add_argument('--plot_dir', required=True)
    parser.add_argument('--min_genes', type=int, default=200)
    parser.add_argument('--min_cells', type=int, default=3)
    parser.add_argument('--max_mito_pct', type=float, default=20.0)
    args = parser.parse_args()

    run_qc(args.input, args.output, args.plot_dir,
           args.min_genes, args.min_cells, args.max_mito_pct)
