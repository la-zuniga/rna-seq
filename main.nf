nextflow.enable.dsl=2

projectDir = "/home/luis/rna-seq"
params.input_dir = "/home/luis/rna-seq/fastqs"
params.reads = "$params.input_dir/*_{1,2}.fastq.gz"
params.outdir = "$projectDir/out"

include { star_index }       from './processes/star_index.nf'
include { fastqc }           from './processes/fastqc.nf'
include { multiqc }          from './processes/multiqc.nf'
include { starsolo }         from './processes/starsolo.nf'
include { merge_counts }     from './processes/merge_counts.nf'
include { scanpy_qc }        from './processes/scanpy_qc.nf'
include { scanpy_preprocess } from './processes/scanpy_preprocess.nf'
include { scanpy_integrate }  from './processes/scanpy_integrate.nf'
include { scanpy_cluster }    from './processes/scanpy_cluster.nf'
include { pseudobulk_de }     from './processes/pseudobulk_de.nf'
include { pathway_analysis }  from './processes/pathway_analysis.nf'

Channel.fromFilePairs(params.reads).set { read_pairs_ch }

workflow {
    genome_ch = Channel.fromPath(params.genome)
    gtf_ch    = Channel.fromPath(params.gtf)
    star_idx  = star_index(genome_ch, gtf_ch)

    fastqcOutput = fastqc(read_pairs_ch)
    multiqcReport = multiqc(fastqcOutput[0].map { sid, qc_dir -> qc_dir }.collect())

    starsoloOutput = starsolo(read_pairs_ch, star_idx)
    merged = merge_counts(starsoloOutput[0].collect())

    qcFiltered = scanpy_qc(merged)
    preprocessed = scanpy_preprocess(qcFiltered[0])
    integrated = scanpy_integrate(preprocessed[0])
    annotated = scanpy_cluster(integrated[0])

    if (params.run_de) {
        deResults = pseudobulk_de(annotated[0])
        pathwayResults = pathway_analysis(deResults[0])
    }
}
