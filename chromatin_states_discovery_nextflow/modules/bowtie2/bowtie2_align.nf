process Bowtie2Alignment
{
    tag "Aligning ${name} with Bowtie2"
    publishDir "$params.outdir/Alignment_logs/",mode:'copy', pattern:'*.log'

    input:
    tuple val(name), val(technology), val(seqmode), path(fastq_files)
    path(genomeIndex)

    output:
    tuple val(name), val(seqmode), val(technology), file("${name}.sam")

    script:
    def cmd = ""
    if (technology == 'ChIP' && seqmode == 'PE') {
        cmd = """
        bowtie2 --sensitive -p 8 -X 2000 -N 0 --trim5 0 --trim3 0 -x ${genomeIndex}/${genomeIndex} \
        -1 ${fastq_files[0]} -2 ${fastq_files[1]} -S ${name}.sam 2> ${name}.bowtie2.log
        """
    } 
    else if (technology == 'ChIP' && seqmode == 'SR') {
        cmd = """
        bowtie2 --sensitive -p 8 -X 2000 -N 0 --trim5 0 --trim3 0 -x ${genomeIndex}/${genomeIndex} \
        -U ${fastq_files} -S ${name}.sam 2> ${name}.bowtie2.log
        """
    }
    else if (technology == 'CNR' && seqmode == 'PE') {
        cmd = """
        bowtie2 --end-to-end --very-sensitive --no-mixed --no-discordant -q \
        --phred33 -I 10 -X 700 -p 8 -x ${genomeIndex}/${genomeIndex} \
        -1 ${fastq_files[0]} -2 ${fastq_files[1]} -S ${name}.sam 2> ${name}.bowtie2.log
        """
    }
    else {
        error "Unsupported alignment configuration: technology=${technology}, seqmode=${seqmode}"
    }

    """
    ${cmd}
    """
}
