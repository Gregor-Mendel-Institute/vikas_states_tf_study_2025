process Sam2Bam
{
    tag "Converting ${name}.sam to sorted BAM"

    input:
    tuple val(name), val(seqmode), val(technology), path(sam_file)

    output:
    tuple val(name), file("${name}.sorted.bam")

    script:
    if(technology == 'ChIP')
    """
    samtools view -S -b ${sam_file} -@ 8 > ${name}.bam
    samtools sort ${name}.bam -@ 8 -o ${name}.sorted.bam
    """

    //size selection step for CNR samples
    else if(technology == 'CNR')
    """
    samtools view -S -b ${sam_file} -@ 8 > ${name}.n.bam
    samtools view -h ${name}.n.bam | awk '\$9 > 150 || \$9 < -150 || \$1 ~ /^@/' | samtools view -bS > ${name}.bam
    samtools sort ${name}.bam -@ 8 -o ${name}.sorted.bam
    """
    else { error "Unsupported technology: ${technology}" }
}
