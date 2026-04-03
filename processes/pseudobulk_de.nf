process pseudobulk_de {
    tag "Pseudobulk DE analysis"
    publishDir "${params.outdir}/pseudobulk_de", mode: 'copy'

    container "${params.containers.scanpy}"

    input:
    path(annotated_h5ad)

    output:
    path("de_results")
    path("de_plots")

    script:
    """
    run_pseudobulk_de.py \
        --input ${annotated_h5ad} \
        --output_dir de_results \
        --plot_dir de_plots \
        --condition_col ${params.condition_col}
    """
}
