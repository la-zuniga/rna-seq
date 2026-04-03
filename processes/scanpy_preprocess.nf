process scanpy_preprocess {
    tag "Scanpy preprocessing"
    publishDir "${params.outdir}/scanpy_preprocess", mode: 'copy'

    container "${params.containers.scanpy}"

    input:
    path(filtered_h5ad)

    output:
    path("preprocessed.h5ad")

    script:
    """
    run_preprocess.py \
        --input ${filtered_h5ad} \
        --output preprocessed.h5ad \
        --n_hvgs ${params.n_hvgs}
    """
}
