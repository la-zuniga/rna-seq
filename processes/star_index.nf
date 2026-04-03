process star_index {
    tag "STAR genome index"
    publishDir "${params.outdir}/star_index", mode: 'copy'

    container "${params.containers.star}"

    input:
    path(genome_fa)
    path(gtf)

    output:
    path("star_genome_index")

    script:
    """
    mkdir star_genome_index

    gunzip -c ${genome_fa} > genome.fa
    gunzip -c ${gtf} > annotations.gtf

    STAR --runMode genomeGenerate \
         --genomeDir star_genome_index \
         --genomeFastaFiles genome.fa \
         --sjdbGTFfile annotations.gtf \
         --runThreadN ${params.star_threads} \
         --limitGenomeGenerateRAM ${params.star_ram}

    rm genome.fa annotations.gtf
    """
}
