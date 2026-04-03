#!/usr/bin/env python3
import argparse
import os
import scanpy as sc
import anndata as ad


def merge_matrices(input_dirs, output_file):
    adatas = []
    for matrix_dir in input_dirs:
        dir_name = os.path.basename(os.path.normpath(matrix_dir))
        parent = os.path.dirname(os.path.normpath(matrix_dir))
        parts = os.path.normpath(matrix_dir).split(os.sep)
        sample_id = None
        for part in parts:
            if '_Solo.out' in part:
                sample_id = part.replace('_Solo.out', '')
                break
        if sample_id is None:
            sample_id = parts[-3] if len(parts) >= 3 else dir_name

        adata = sc.read_10x_mtx(matrix_dir, var_names='gene_symbols', cache=False)
        adata.obs['sample_id'] = sample_id
        adata.var_names_make_unique()
        adatas.append(adata)
        print(f"  {sample_id}: {adata.n_obs} cells x {adata.n_vars} genes")

    combined = ad.concat(
        adatas,
        join='outer',
        label='sample_id',
        keys=[a.obs['sample_id'].iloc[0] for a in adatas]
    )
    combined.obs_names_make_unique()

    print(f"\nMerged: {combined.n_obs} total cells x {combined.n_vars} genes")
    print(f"Samples: {combined.obs['sample_id'].nunique()}")

    combined.write_h5ad(output_file)
    print(f"Written to {output_file}")


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--input_dirs', nargs='+', required=True)
    parser.add_argument('--output', required=True)
    args = parser.parse_args()

    merge_matrices(args.input_dirs, args.output)
