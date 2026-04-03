process scanpy_cluster {
    tag "Scanpy clustering"
    publishDir "${params.outdir}/scanpy_cluster", mode: 'copy'

    container "${params.containers.scanpy}"

    input:
    path(integrated_h5ad)

    output:
    path("annotated.h5ad")
    path("cluster_plots")
    path("marker_genes")

    script:
    """
    run_clustering.py \
        --input ${integrated_h5ad} \
        --output annotated.h5ad \
        --plot_dir cluster_plots \
        --marker_dir marker_genes \
        --resolution ${params.leiden_res}
    """
}
