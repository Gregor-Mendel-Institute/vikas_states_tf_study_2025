process overlapenrichment
{
    tag "Calculating enrichment for ${name}"
    publishDir "${params.outdir}/chromhmm/overlapenrichment", mode: 'copy'

    input:
    tuple val(name), path(segments_file), path(feature_dir) 

    output:
    path "${name}_overlap.*" 

    script:
    """
    java -Xmx16G -jar /ChromHMM/ChromHMM.jar OverlapEnrichment \
    -color 0,0,0 "${segments_file}" "${feature_dir}" "${name}_overlap"
    """
}
