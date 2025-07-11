process create_bowtie2_index {
    tag "Creating index for custom genome"
    publishDir "${params.outdir}/genome_index", mode: 'copy'

    input:
    path reference_fasta
    val index_basename

    output:
    path "${index_basename}"

    when:
    params.reference_fasta != "null"

    script:
    """
    mkdir -p ${index_basename}
    bowtie2-build --threads ${task.cpus} ${reference_fasta} ${index_basename}/${index_basename}
    """
}
