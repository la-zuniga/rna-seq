process fastqc {
    tag "FastQC on $sample_id"
    publishDir "${params.outdir}/${sample_id}/fastqc", mode: 'copy'

    container "${params.containers.fastqc}"

    input:
    tuple val(sample_id), path(reads)

    output:
    tuple val(sample_id), path("fastqc_results")

    script:
    """
    mkdir fastqc_results
    fastqc ${reads} --outdir fastqc_results
    """
}
