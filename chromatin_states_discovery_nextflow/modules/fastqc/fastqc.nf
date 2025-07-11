process FASTQC {
    tag "Running FASTQC on ${name}"
    publishDir "${params.outdir}/First_Fastqc_Reports", mode: 'copy'

    input:
    tuple val(name), val(technology), val(seqmode), path(fastq_file)

    output:
    path "${name}*fastqc.{zip,html}"

    script:
    """
    fastqc ${fastq_file}
    """
}
