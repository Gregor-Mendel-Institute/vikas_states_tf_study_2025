def handleFastqInputs() {
    if (params.fastq_file_type == "Only_SR") {
        fq_sr = Channel.fromPath(params.fqs_sr, checkIfExists: true)
        fq_sr_sample_Info = Channel.fromPath(params.fq_samples).splitCsv(sep: "\t").map { row -> tuple(row[0], row[1], row[2]) }
        return fq_sr_sample_Info.join(fq_sr)
    } else if (params.fastq_file_type == "Only_PE") {
        fq_pe = Channel.fromFilePairs(params.fqs_pe, checkIfExists: true)
        fq_pe_sample_Info = Channel.fromPath(params.fq_samples).splitCsv(sep: "\t").map { row -> tuple(row[0], row[1], row[2]) }
        return fq_pe_sample_Info.join(fq_pe)
    } else if (params.fastq_file_type == "SR_and_PE") {
        // Single-end processing
        fq_sr = Channel.fromPath(params.fqs_sr, checkIfExists: true)
        fq_sr_sample_Info = Channel.fromPath(params.fq_samples).splitCsv(sep: "\t").map { row -> tuple(row[0], row[1], row[2]) }
        fq_comb_sr_Info = fq_sr_sample_Info.join(fq_sr)

        // Paired-end processing
        fq_pe = Channel.fromFilePairs(params.fqs_pe, checkIfExists: true)
        fq_pe_sample_Info = Channel.fromPath(params.fq_samples).splitCsv(sep: "\t").map { row -> tuple(row[0], row[1], row[2]) }
        fq_comb_pe_Info = fq_pe_sample_Info.join(fq_pe)

        // Combine both
        return fq_comb_sr_Info.mix(fq_comb_pe_Info)
    } else {
        error "Invalid fastq_file_type parameter. Must be one of: Only_SR, Only_PE, SR_and_PE"
    }
}

def handleBamInputs() {
    bam_files = Channel.fromPath(params.bams, checkIfExists: true)
                       .map { file -> tuple(file.baseName, file) }
    bam_sample_Info = Channel.fromPath(params.bam_samples)
                             .splitCsv(sep: "\t")
                             .map { row -> tuple(row[0], row[1], row[2]) }
    return bam_sample_Info.join(bam_files)
}
