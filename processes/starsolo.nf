process starsolo {
    tag "STARsolo on $sample_id"
    publishDir "${params.outdir}/${sample_id}/starsolo", mode: 'copy'

    container "${params.containers.star}"

    input:
    tuple val(sample_id), path(reads)
    path(star_index)

    output:
    tuple val(sample_id), path("${sample_id}_Solo.out/Gene/filtered")
    path("${sample_id}_Aligned.sortedByCoord.out.bam")
    path("${sample_id}_Log.final.out")

    script:
    """
    STAR --genomeDir ${star_index} \
         --readFilesIn ${reads[1]} ${reads[0]} \
         --readFilesCommand zcat \
         --soloType ${params.solo_type} \
         --soloCBwhitelist ${params.whitelist} \
         --soloCBlen ${params.solo_cb_len} \
         --soloUMIlen ${params.solo_umi_len} \
         --outSAMtype BAM SortedByCoordinate \
         --outSAMattributes NH HI nM AS CR UR CB UB GX GN sS sQ sM \
         --soloFeatures Gene \
         --outFileNamePrefix ${sample_id}_ \
         --runThreadN ${params.star_threads}
    """
}
