#!/usr/bin/env python3
import argparse
import os
import pandas as pd
import scanpy as sc
import celltypist
from celltypist import models
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt


def run_clustering(input_file, output_file, plot_dir, marker_dir, resolution):
    os.makedirs(plot_dir, exist_ok=True)
    os.makedirs(marker_dir, exist_ok=True)

    adata = sc.read_h5ad(input_file)
    print(f"Loaded: {adata.n_obs} cells x {adata.n_vars} genes")

    sc.tl.leiden(adata, resolution=resolution, flavor='igraph')
    n_clusters = adata.obs['leiden'].nunique()
    print(f"Leiden (res={resolution}): {n_clusters} clusters")

    sc.tl.rank_genes_groups(adata, 'leiden', method='wilcoxon')

    result = adata.uns['rank_genes_groups']
    for cluster in adata.obs['leiden'].unique():
        markers_df = pd.DataFrame({
            'gene': result['names'][cluster],
            'logfoldchange': result['logfoldchanges'][cluster],
            'pval': result['pvals'][cluster],
            'pval_adj': result['pvals_adj'][cluster],
            'score': result['scores'][cluster]
        })
        markers_df.to_csv(
            os.path.join(marker_dir, f'cluster_{cluster}_markers.tsv'),
            sep='\t', index=False
        )

    model = models.download_model('Immune_All_Low.pkl')
    model = celltypist.models.Model.load(model='Immune_All_Low.pkl')

    if adata.raw is not None:
        adata_for_ct = adata.raw.to_adata()
    else:
        adata_for_ct = adata.copy()

    predictions = celltypist.annotate(
        adata_for_ct,
        model=model,
        majority_voting=True
    )
    adata.obs['cell_type'] = predictions.predicted_labels['majority_voting']

    for ct, count in adata.obs['cell_type'].value_counts().items():
        print(f"  {ct}: {count} ({100 * count / adata.n_obs:.1f}%)")

    sc.pl.umap(adata, color='leiden', show=False, legend_loc='on data')
    plt.savefig(os.path.join(plot_dir, 'umap_leiden.png'), dpi=150, bbox_inches='tight')
    plt.close()

    sc.pl.umap(adata, color='cell_type', show=False)
    plt.savefig(os.path.join(plot_dir, 'umap_cell_type.png'), dpi=150, bbox_inches='tight')
    plt.close()

    sc.pl.rank_genes_groups_dotplot(adata, n_genes=5, show=False)
    plt.savefig(os.path.join(plot_dir, 'dotplot_markers.png'), dpi=150, bbox_inches='tight')
    plt.close()

    adata.write_h5ad(output_file)
    print(f"\nWritten to {output_file}")


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--input', required=True)
    parser.add_argument('--output', required=True)
    parser.add_argument('--plot_dir', required=True)
    parser.add_argument('--marker_dir', required=True)
    parser.add_argument('--resolution', type=float, default=1.0)
    args = parser.parse_args()

    run_clustering(args.input, args.output, args.plot_dir,
                   args.marker_dir, args.resolution)
