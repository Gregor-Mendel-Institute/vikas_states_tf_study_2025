process TrimAdaptors {
    tag "Running TrimGalore! on ${name}"
    publishDir "${params.outdir}/FilterFastq", mode: 'copy', pattern: '*.{html,zip,txt}'

    input:
    tuple val(name), val(technology), val(seqmode), path(fastq_file)

    output:
    tuple val(name), val(technology), val(seqmode), path("${name}_trimmed*fastq")

    script:
    if (seqmode == 'SR') {
        """
        trim_galore --dont_gzip --stringency 1 \
        --fastqc --length 5 ${fastq_file} --basename ${name}

        mv ${name}_trimmed.fq ${name}_trimmed.fastq
        """
    } 
    else if (seqmode == 'PE') {
        """
        trim_galore --dont_gzip --stringency 1 \
        --fastqc --length 5 --paired ${fastq_file[0]} ${fastq_file[1]} \
        --basename ${name}

        mv "${name}_val_1.fq" "${name}_trimmed_1.fastq"
        mv "${name}_val_2.fq" "${name}_trimmed_2.fastq"
        """
    } else {
        error "Unknown sequencing mode: ${seqmode}. Must be 'SR' or 'PE'"
    }
}
