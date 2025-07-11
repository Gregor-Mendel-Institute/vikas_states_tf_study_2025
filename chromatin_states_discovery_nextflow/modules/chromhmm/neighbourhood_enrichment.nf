process NeighbourhoodEnrichment
{
    tag "Running NE-analysis for ${name}"
    publishDir "${params.outdir}/chromhmm/Neighbourhood_Enrichment/TSS", mode: 'copy', pattern:'*_TSS.*'
    publishDir "${params.outdir}/chromhmm/Neighbourhood_Enrichment/TTS", mode: 'copy', pattern:'*_TTS.*'

    input:
    tuple val(name), path(segments_bed), path(TSS), file(TTS) 

    output:
    path "${name}_TSS.*"
    path "${name}_TTS.*"

    script:
    """
    java -Xmx16G -jar /ChromHMM/ChromHMM.jar NeighborhoodEnrichment \
    -t "Position relative to TSS" "${segments_bed}" "${TSS}" "${name}_TSS"

    java -Xmx16G -jar /ChromHMM/ChromHMM.jar NeighborhoodEnrichment \
    -t "Position relative to TTS" "${segments_bed}" "${TTS}" "${name}_TTS"
    """
}
