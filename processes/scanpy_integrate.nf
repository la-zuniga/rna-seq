process scanpy_integrate {
    tag "Scanpy integration"
    publishDir "${params.outdir}/scanpy_integrate", mode: 'copy'

    container "${params.containers.scanpy}"

    input:
    path(preprocessed_h5ad)

    output:
    path("integrated.h5ad")
    path("integration_plots")

    script:
    """
    run_integration.py \
        --input ${preprocessed_h5ad} \
        --output integrated.h5ad \
        --plot_dir integration_plots \
        --n_pcs ${params.n_pcs} \
        --batch_key sample_id
    """
}
