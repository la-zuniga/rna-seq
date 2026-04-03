#!/usr/bin/env python3
import argparse
import scanpy as sc


def run_preprocess(input_file, output_file, n_hvgs):
    adata = sc.read_h5ad(input_file)
    print(f"Loaded: {adata.n_obs} cells x {adata.n_vars} genes")

    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    adata.raw = adata.copy()

    sc.pp.highly_variable_genes(
        adata,
        n_top_genes=n_hvgs,
        batch_key='sample_id',
        flavor='seurat_v3',
        subset=False
    )
    n_hvg = adata.var['highly_variable'].sum()
    print(f"{n_hvg} HVGs selected")

    adata = adata[:, adata.var['highly_variable']].copy()
    sc.pp.scale(adata, max_value=10)

    adata.write_h5ad(output_file)
    print(f"Output: {adata.n_obs} cells x {adata.n_vars} genes -> {output_file}")


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--input', required=True)
    parser.add_argument('--output', required=True)
    parser.add_argument('--n_hvgs', type=int, default=2000)
    args = parser.parse_args()

    run_preprocess(args.input, args.output, args.n_hvgs)
