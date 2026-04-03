process pathway_analysis {
    tag "Pathway analysis"
    publishDir "${params.outdir}/pathway_analysis", mode: 'copy'

    container "${params.containers.scanpy}"

    input:
    path(de_results)

    output:
    path("pathway_results")
    path("pathway_plots")

    script:
    """
    run_pathway_analysis.py \
        --input_dir ${de_results} \
        --output_dir pathway_results \
        --plot_dir pathway_plots
    """
}
