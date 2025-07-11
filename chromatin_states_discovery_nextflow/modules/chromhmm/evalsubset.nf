process EvalSubset
{
    tag "Evaluating subset for ${order_name} of ${num}-state model"
    publishDir "${params.outdir}/chromhmm/evalSubset/${order_name}/", mode:'copy', pattern:"*_${order_name}_*"

    input:
    tuple val(num), path(model), path(segments_dir), path(binarized_chrs), val(order_name), val(order_string)

    output:
    path "${num}_${order_name}_evasub.*"

    script:
    """
    java -Xmx16G -jar /ChromHMM/ChromHMM.jar EvalSubset \
    ${model} ${binarized_chrs} ${segments_dir} "${num}_${order_name}_evasub" \
    ${order_string}
    """
}
