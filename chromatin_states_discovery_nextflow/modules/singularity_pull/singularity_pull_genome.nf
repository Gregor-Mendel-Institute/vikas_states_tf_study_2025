process pullGenomeIndex {
    tag "Pulling pre-built index for ${params.reference_genome}"
    publishDir "${params.outdir}/genome_index", mode: 'copy'

    input:
    val index

    output:
    path "${params.reference_genome}/*"

    when:
    params.reference_fasta == "null" && params.reference_genome != "null"

    script:
    """
    mkdir -p ${params.reference_genome}
    singularity run ${index}
    """
}
