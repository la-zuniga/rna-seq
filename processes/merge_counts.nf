process merge_counts {
    tag "Merge count matrices"
    publishDir "${params.outdir}/merged", mode: 'copy'

    container "${params.containers.scanpy}"

    input:
    path(sample_dirs)

    output:
    path("merged_counts.h5ad")

    script:
    """
    merge_to_h5ad.py \
        --input_dirs ${sample_dirs} \
        --output merged_counts.h5ad
    """
}
