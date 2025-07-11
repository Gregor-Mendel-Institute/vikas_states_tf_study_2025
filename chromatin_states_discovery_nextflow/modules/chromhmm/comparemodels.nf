process comparemodels
{
    tag "Comparemodel for ${max_numstates}-state model"
    publishDir "${params.outdir}/chromhmm/comparemodels", mode: 'copy'

    input:
    path emission_files
    val max_numstates

    output:
    path ("${max_numstates}_comparemodels.*") 

    script:
    """
    java -Xmx16G -jar /ChromHMM/ChromHMM.jar CompareModels -color 0,255,0 \
    emissions_${max_numstates}.txt . ${max_numstates}_comparemodels
    """
}

