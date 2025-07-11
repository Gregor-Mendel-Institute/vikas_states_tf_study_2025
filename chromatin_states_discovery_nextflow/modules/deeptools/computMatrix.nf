process computeMatrix
{
    tag "computing matrix for ${mark}"
    publishDir "$params.outdir/deeptools_matrix/", mode:'copy'

    input:
    tuple val(mark), path(bigwig_files), path(regions_file)

    output:
    tuple val(mark), path("${mark}_matrix.gz")

    script:
    """
    computeMatrix reference-point --referencePoint TSS \
    -a 3000 -b 1000 -R ${regions_file} \
    -S ${bigwig_files.join(' ')} --skipZeros --missingDataAsZero \
    --outFileName ${mark}_matrix.gz
    """
}
