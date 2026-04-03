#!/usr/bin/env python3
import argparse
import os
import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt


def aggregate_pseudobulk(adata, cell_type_col, sample_col, condition_col):
    pseudobulk_data = {}

    for cell_type in adata.obs[cell_type_col].unique():
        ct_mask = adata.obs[cell_type_col] == cell_type
        ct_adata = adata[ct_mask]

        counts_list = []
        meta_list = []

        for sample in ct_adata.obs[sample_col].unique():
            sample_mask = ct_adata.obs[sample_col] == sample
            sample_adata = ct_adata[sample_mask]

            if sample_adata.n_obs < 10:
                continue

            if adata.raw is not None:
                raw_subset = adata.raw[ct_mask & (adata.obs[sample_col] == sample)]
                counts = np.array(raw_subset.X.sum(axis=0)).flatten()
                gene_names = adata.raw.var_names
            else:
                counts = np.array(sample_adata.X.sum(axis=0)).flatten()
                gene_names = adata.var_names

            counts_list.append(counts)
            condition = sample_adata.obs[condition_col].iloc[0]
            meta_list.append({'sample': sample, condition_col: condition})

        if len(counts_list) < 4:
            print(f"  {cell_type}: skipped ({len(counts_list)} samples, need >= 4)")
            continue

        count_df = pd.DataFrame(
            counts_list,
            index=[m['sample'] for m in meta_list],
            columns=gene_names
        ).astype(int)
        meta_df = pd.DataFrame(meta_list).set_index('sample')

        pseudobulk_data[cell_type] = (count_df, meta_df)
        print(f"  {cell_type}: {len(counts_list)} samples")

    return pseudobulk_data


def run_deseq2(count_df, meta_df, condition_col):
    from pydeseq2.dds import DeseqDataSet
    from pydeseq2.ds import DeseqStats

    keep = count_df.sum(axis=0) >= 10
    count_df = count_df.loc[:, keep]

    dds = DeseqDataSet(
        counts=count_df,
        metadata=meta_df,
        design_factors=condition_col
    )
    dds.deseq2()

    conditions = meta_df[condition_col].unique()
    stat_res = DeseqStats(
        dds,
        contrast=[condition_col, conditions[0], conditions[1]]
    )
    stat_res.summary()

    return stat_res.results_df


def volcano_plot(results_df, cell_type, output_path, pval_thresh=0.05, lfc_thresh=1.0):
    results_df = results_df.dropna(subset=['padj', 'log2FoldChange'])

    fig, ax = plt.subplots(figsize=(8, 6))
    colors = np.where(
        (results_df['padj'] < pval_thresh) & (abs(results_df['log2FoldChange']) > lfc_thresh),
        'red', 'grey'
    )
    ax.scatter(
        results_df['log2FoldChange'],
        -np.log10(results_df['padj']),
        c=colors, alpha=0.5, s=5
    )
    ax.axhline(-np.log10(pval_thresh), color='black', linestyle='--', linewidth=0.5)
    ax.axvline(-lfc_thresh, color='black', linestyle='--', linewidth=0.5)
    ax.axvline(lfc_thresh, color='black', linestyle='--', linewidth=0.5)
    ax.set_xlabel('log2 Fold Change')
    ax.set_ylabel('-log10 adjusted p-value')
    ax.set_title(f'Volcano Plot: {cell_type}')
    plt.tight_layout()
    plt.savefig(output_path, dpi=150)
    plt.close()


def run_pseudobulk_de(input_file, output_dir, plot_dir, condition_col):
    os.makedirs(output_dir, exist_ok=True)
    os.makedirs(plot_dir, exist_ok=True)

    adata = sc.read_h5ad(input_file)
    print(f"Loaded: {adata.n_obs} cells x {adata.n_vars} genes")

    pseudobulk_data = aggregate_pseudobulk(
        adata, 'cell_type', 'sample_id', condition_col
    )

    for cell_type, (count_df, meta_df) in pseudobulk_data.items():
        safe_name = cell_type.replace(' ', '_').replace('/', '_')
        print(f"\n  {cell_type}:")

        try:
            results_df = run_deseq2(count_df, meta_df, condition_col)
            results_df.to_csv(
                os.path.join(output_dir, f'{safe_name}_de_results.tsv'),
                sep='\t'
            )

            n_sig = (results_df['padj'] < 0.05).sum()
            print(f"    {len(results_df)} genes tested, {n_sig} significant (padj < 0.05)")

            volcano_plot(
                results_df, cell_type,
                os.path.join(plot_dir, f'{safe_name}_volcano.png')
            )

        except Exception as e:
            print(f"    Error: {e}")

    print(f"\nResults written to {output_dir}")


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--input', required=True)
    parser.add_argument('--output_dir', required=True)
    parser.add_argument('--plot_dir', required=True)
    parser.add_argument('--condition_col', default='condition')
    args = parser.parse_args()

    run_pseudobulk_de(args.input, args.output_dir, args.plot_dir, args.condition_col)
