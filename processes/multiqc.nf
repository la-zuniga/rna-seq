process multiqc {
    tag "MultiQC aggregation"
    publishDir "${params.outdir}/multiqc", mode: 'copy'

    container "${params.containers.multiqc}"

    input:
    path(qc_dirs)

    output:
    path("MultiQC_report.html")

    script:
    """
    multiqc ${qc_dirs} -o . -n MultiQC_report.html
    """
}
