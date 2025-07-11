def validateGenomeSelection() {
    if (params.reference_fasta == "null" && params.reference_genome == "null") {
        error "Neither reference FASTA nor reference genome specified. Please provide one."
    }

    if (params.reference_genome != "null" && !params.GENOME_INDICES.containsKey(params.reference_genome)) {
        error """Invalid reference genome: ${params.reference_genome}
                Available genomes: ${params.GENOME_INDICES.keySet().join(', ')}"""
    }

    log.info "Validation successful: Reference genome or FASTA provided."
}
