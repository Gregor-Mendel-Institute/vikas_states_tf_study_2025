process Bam2Fastq
{
    tag "BAM-to-FASTQ conversion for ${name}"

    input:
    tuple val(name), val(technology), val(seqmode), path(bam)

    output:
    tuple val(name), val(technology), val(seqmode), path("${name}*.fastq")

    script:
    if (seqmode == 'SR')
    """
    samtools sort -n -T ${name}.tmp ${bam} | bedtools bamtofastq -i stdin -fq ${name}.fastq
    """
    else if (seqmode == 'PE')
    """
    samtools sort -n -T ${name}.tmp ${bam} | bedtools bamtofastq -i stdin -fq ${name}_1.fastq -fq2 ${name}_2.fastq
    """
    else { error "Unknown sequencing mode: ${seqmode}. Must be 'SR' or 'PE'" }
}
