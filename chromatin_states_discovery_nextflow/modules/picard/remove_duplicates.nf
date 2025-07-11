process rmdup
{
    tag "Removing duplicates from ${name}.sorted.bam"
    publishDir "$params.outdir/Aligned_dedup_BAMS/", mode:'copy', pattern:'*.dedup.metrics.txt'

    input:
    tuple val(name), path(sorted_bam)

    output:
    tuple val(name), file("${name}.dedup.bam")

    script:
    """
    java -jar /usr/picard/picard.jar MarkDuplicates \
    I=${sorted_bam} \
    O=${name}.dedup.bam \
    M=${name}.dedup.metrics.txt \
    ASSUME_SORTED=true \
    REMOVE_DUPLICATES=true
    """
}
