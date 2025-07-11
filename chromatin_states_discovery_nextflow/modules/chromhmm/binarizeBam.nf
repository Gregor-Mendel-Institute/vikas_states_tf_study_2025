process binarize
{
    tag "Binarizing BAMs"
    publishDir "$params.outdir/chromhmm/", mode:'copy'

    input:
    path (BAMs_BAIs_dir)
    path (cellmarkfiletable)
    path (chrlengthfile)

    output:
    path (binarized_files)

    script:
    """
    java -Xmx16G -jar /ChromHMM/ChromHMM.jar BinarizeBam \
    -gzip -mixed \
    ${chrlengthfile} ${BAMs_BAIs_dir} ${cellmarkfiletable} "binarized_files/"
    """
}
