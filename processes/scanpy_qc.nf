process scanpy_qc {
    tag "Scanpy QC filtering"
    publishDir "${params.outdir}/scanpy_qc", mode: 'copy'

    container "${params.containers.scanpy}"

    input:
    path(merged_h5ad)

    output:
    path("qc_filtered.h5ad")
    path("qc_plots")

    script:
    """
    run_qc.py \
        --input ${merged_h5ad} \
        --output qc_filtered.h5ad \
        --plot_dir qc_plots \
        --min_genes ${params.min_genes} \
        --min_cells ${params.min_cells} \
        --max_mito_pct ${params.max_mito_pct}
    """
}
