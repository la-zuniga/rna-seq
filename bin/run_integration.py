#!/usr/bin/env python3
import argparse
import os
import scanpy as sc
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt


def run_integration(input_file, output_file, plot_dir, n_pcs, batch_key):
    os.makedirs(plot_dir, exist_ok=True)

    adata = sc.read_h5ad(input_file)
    print(f"Loaded: {adata.n_obs} cells x {adata.n_vars} genes")

    sc.tl.pca(adata, n_comps=n_pcs)

    n_samples = adata.obs[batch_key].nunique()
    if n_samples > 1:
        sc.external.pp.harmony_integrate(adata, key=batch_key)
        use_rep = 'X_pca_harmony'
        print(f"Harmony: {n_samples} batches on '{batch_key}'")
    else:
        use_rep = 'X_pca'
        print("Single sample, skipping Harmony")

    sc.pp.neighbors(adata, use_rep=use_rep)
    sc.tl.umap(adata)

    sc.pl.umap(adata, color=batch_key, show=False)
    plt.savefig(os.path.join(plot_dir, 'umap_by_sample.png'), dpi=150, bbox_inches='tight')
    plt.close()

    sc.pl.umap(adata, color='n_genes_by_counts', show=False)
    plt.savefig(os.path.join(plot_dir, 'umap_by_ngenes.png'), dpi=150, bbox_inches='tight')
    plt.close()

    sc.pl.umap(adata, color='total_counts', show=False)
    plt.savefig(os.path.join(plot_dir, 'umap_by_counts.png'), dpi=150, bbox_inches='tight')
    plt.close()

    adata.write_h5ad(output_file)
    print(f"Written to {output_file}")


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--input', required=True)
    parser.add_argument('--output', required=True)
    parser.add_argument('--plot_dir', required=True)
    parser.add_argument('--n_pcs', type=int, default=50)
    parser.add_argument('--batch_key', default='sample_id')
    args = parser.parse_args()

    run_integration(args.input, args.output, args.plot_dir,
                    args.n_pcs, args.batch_key)
