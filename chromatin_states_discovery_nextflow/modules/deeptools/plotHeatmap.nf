process plotHeatmap
{
    tag "Plotting heatmap for ${mark}"
    publishDir "${params.outdir}/Heatmaps/", mode:'copy'

    input:
    tuple val(mark), path(matrix)

    output:
    path ("${mark}_heatmap.pdf")

    script:
    """
    plotHeatmap -m ${matrix} -out ${mark}_heatmap.pdf \
    --colorMap binary --labelRotation 45 --kmeans 4
    """
}
