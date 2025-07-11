process index_rmdup
{
    tag "Indexing the deduplicated file ${name}.sorted.bam"
    publishDir "$params.outdir/Aligned_dedup_BAMS/", mode:'copy'

    input:
    tuple val(name), path(dedup_bam)

    output:
    tuple val(name), path("${name}.dedup.sorted.bam"), path("${name}.dedup.sorted.bam.bai")

    script:
    """
    samtools sort ${dedup_bam} -@ 16 > ${name}.dedup.sorted.bam
    samtools index ${name}.dedup.sorted.bam -@ 16
    """
}
