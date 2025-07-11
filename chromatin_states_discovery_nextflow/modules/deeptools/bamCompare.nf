process bamCompare
{
    tag "${test} vs ${control} -> ${out}.bw"
    publishDir path: {
        if (out.contains('inputNormalized')) {
            return "$params.outdir/BigWigs/Input-Normalized/"
        } else if (out.contains('H3Normalized')) {
            return "$params.outdir/BigWigs/H3-Normalized/"
        } else {
            return "$params.outdir/BigWigs/"
        }
    }, mode: 'copy'

    input:
    tuple val(test), path(test_bam), path(test_bai), val(control), path(control_bam), path(control_bai), val(out)

    output:
    tuple val(out), path("${out}.bw")

    script:
    """
    bamCompare -b1 ${test_bam} -b2 ${control_bam} -o ${out}.bw -p max
    """
}
