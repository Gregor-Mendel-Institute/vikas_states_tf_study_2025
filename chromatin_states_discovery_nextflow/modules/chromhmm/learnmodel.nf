process learnmodel
{
    tag "ChromHMM Learning model for ${numstates}"
    publishDir "$params.outdir/chromhmm/", mode:'copy', pattern:'learnmodel/*'

    input:
    tuple path (binarized_dir), val (assembly_name), val (numstates)

    output:
    path ("learnmodel/emissions_${numstates}.txt"), emit: emissions
    path ("learnmodel/*${numstates}_segments.bed.gz"), emit: segments
    tuple val("${numstates}"), path("learnmodel/model_${numstates}.txt"), emit: allmodels

    script:
    """
    java -Xmx16G -jar /ChromHMM/ChromHMM.jar LearnModel -p 0 \
    -gzip -color 0,255,0 ${binarized_dir} "learnmodel/" \
    ${numstates} ${assembly_name}
    """
}
