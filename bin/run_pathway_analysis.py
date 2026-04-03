#!/usr/bin/env python3
import argparse
import os
import glob
import numpy as np
import pandas as pd
import gseapy
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt


def run_gsea_for_celltype(de_file, cell_type, output_dir, plot_dir):
    de_df = pd.read_csv(de_file, sep='\t', index_col=0)
    de_df = de_df.dropna(subset=['pvalue', 'log2FoldChange'])

    de_df['rank_metric'] = -np.log10(de_df['pvalue'] + 1e-300) * np.sign(de_df['log2FoldChange'])
    ranked = de_df['rank_metric'].sort_values(ascending=False)

    safe_name = cell_type.replace(' ', '_').replace('/', '_')
    ct_output = os.path.join(output_dir, safe_name)
    os.makedirs(ct_output, exist_ok=True)

    try:
        result = gseapy.prerank(
            rnk=ranked,
            gene_sets='MSigDB_Hallmark_2020',
            outdir=ct_output,
            min_size=15,
            max_size=500,
            permutation_num=1000,
            seed=42,
            verbose=False
        )

        if result.res2d is not None and len(result.res2d) > 0:
            result.res2d.to_csv(
                os.path.join(ct_output, f'{safe_name}_gsea_results.tsv'),
                sep='\t'
            )
            n_sig = (result.res2d['FDR q-val'] < 0.25).sum()
            print(f"  {cell_type}: {len(result.res2d)} pathways, {n_sig} significant (FDR < 0.25)")
            return result.res2d
        else:
            print(f"  {cell_type}: no enrichment results")
            return None

    except Exception as e:
        print(f"  {cell_type}: GSEA error -- {e}")
        return None


def summary_heatmap(all_results, plot_dir):
    if not all_results:
        return

    nes_data = {}
    for cell_type, res_df in all_results.items():
        if res_df is not None:
            for _, row in res_df.iterrows():
                pathway = row['Term']
                if pathway not in nes_data:
                    nes_data[pathway] = {}
                nes_data[pathway][cell_type] = row['NES']

    if not nes_data:
        return

    heatmap_df = pd.DataFrame(nes_data).T.fillna(0)

    top_pathways = heatmap_df.abs().max(axis=1).nlargest(20).index
    heatmap_df = heatmap_df.loc[top_pathways]

    fig, ax = plt.subplots(figsize=(max(8, len(heatmap_df.columns) * 1.5), 10))
    im = ax.imshow(heatmap_df.values, cmap='RdBu_r', aspect='auto', vmin=-3, vmax=3)
    ax.set_xticks(range(len(heatmap_df.columns)))
    ax.set_xticklabels(heatmap_df.columns, rotation=45, ha='right')
    ax.set_yticks(range(len(heatmap_df.index)))
    ax.set_yticklabels(heatmap_df.index, fontsize=8)
    plt.colorbar(im, ax=ax, label='NES')
    ax.set_title('Pathway Enrichment (NES) Across Cell Types')
    plt.tight_layout()
    plt.savefig(os.path.join(plot_dir, 'pathway_heatmap.png'), dpi=150)
    plt.close()
    print(f"Heatmap saved to {plot_dir}/pathway_heatmap.png")


def run_pathway_analysis(input_dir, output_dir, plot_dir):
    os.makedirs(output_dir, exist_ok=True)
    os.makedirs(plot_dir, exist_ok=True)

    de_files = sorted(glob.glob(os.path.join(input_dir, '*_de_results.tsv')))
    print(f"Found {len(de_files)} DE result files")

    all_results = {}
    for de_file in de_files:
        cell_type = os.path.basename(de_file).replace('_de_results.tsv', '')
        result = run_gsea_for_celltype(de_file, cell_type, output_dir, plot_dir)
        all_results[cell_type] = result

    summary_heatmap(all_results, plot_dir)
    print(f"Pathway results in {output_dir}")


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--input_dir', required=True)
    parser.add_argument('--output_dir', required=True)
    parser.add_argument('--plot_dir', required=True)
    args = parser.parse_args()

    run_pathway_analysis(args.input_dir, args.output_dir, args.plot_dir)
