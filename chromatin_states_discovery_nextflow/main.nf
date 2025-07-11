#!/usr/bin/env nextflow

nextflow.enable.dsl=2

// Import modules and functions

include { validateGenomeSelection } from './modules/functions/genome_validation'
include { handleFastqInputs       } from './modules/functions/file_handler_functions'
include { handleBamInputs         } from './modules/functions/file_handler_functions'

include { create_bowtie2_index    } from './modules/bowtie2/bowtie2_built'
include { pullGenomeIndex         } from './modules/singularity_pull/singularity_pull_genome'
include { Bam2Fastq               } from './modules/samtools/bam2fastq'
include { FASTQC                  } from './modules/fastqc/fastqc'
include { TrimAdaptors            } from './modules/trimgalore/trim_galore'
include { Bowtie2Alignment        } from './modules/bowtie2/bowtie2_align'
include { Sam2Bam                 } from './modules/samtools/sam2bam'
include { rmdup                   } from './modules/picard/remove_duplicates'
include { index_rmdup             } from './modules/samtools/sort_index'
include { binarize                } from './modules/chromhmm/binarizeBam'
include { learnmodel              } from './modules/chromhmm/learnmodel'
include { comparemodels           } from './modules/chromhmm/comparemodels'
include { overlapenrichment       } from './modules/chromhmm/overlapenrichment'
include { NeighbourhoodEnrichment } from './modules/chromhmm/neighbourhood_enrichment'
include { EvalSubset              } from './modules/chromhmm/evalsubset'
include { bamCompare              } from './modules/deeptools/bamCompare'
include { computeMatrix           } from './modules/deeptools/computMatrix'
include { plotHeatmap             } from './modules/deeptools/plotHeatmap'

workflow {
            // Print pipeline information
            log.info"""
            
===============================================================================================================

                  C h I P - S E Q   A L I G N M E N T   N F -  P I P E L I N E

Author: Vikas Shukla <vikas.shukla@gmi.oeaw.ac.at>

===============================================================================================================

        File Types                      : ${params.file_types}
        sample Information (fastq)      : ${params.fq_samples}
        sample Information (bam)        : ${params.bam_samples}
        Input (SR fastq)                : ${params.fqs_sr}
        Input (PE fastq)                : ${params.fqs_pe}
        Input (BAM FILES)               : ${params.bams}
        Reference Genome                : ${params.reference_genome}
        output Directory                : ${params.outdir}
            """

            // Genome validation
            validateGenomeSelection()

            // Create channels for custom reference data
            if (params.reference_fasta != "null") {
                reference_fasta_ch = Channel.fromPath(params.reference_fasta, checkIfExists: true)
                // Assembly name is required for learnmodel
                assembly_name = Channel.value(params.index_basename)
                index_basename_ch = Channel.value(params.index_basename)

                // Create Bowtie2 index for provided custom reference
                bw2index = create_bowtie2_index(reference_fasta_ch, index_basename_ch)
            }

            // Pull reference genome from containers for existing genomes
            else if (params.reference_fasta == "null" && params.reference_genome != "null") {
                index_ch = Channel.value(params.GENOME_INDICES[params.reference_genome].index)
                // Assembly name is required for learnmodel
                assembly_name = Channel.value(params.reference_genome)
                bw2index = pullGenomeIndex(index_ch)
            }

            else {
                error "No reference genome or FASTA provided. Cannot proceed!"
            }

            // Handle input files based on input file_types
            def input_fastq_ch = null
            def input_bam_ch = null

            if (params.file_types == "fastq") {
                input_fastq_ch = handleFastqInputs()
            } 
            else if (params.file_types == "bam") {
                input_bam_ch = handleBamInputs()
            } 
            else if (params.file_types == "bam_and_fastq") {
                input_fastq_ch = handleFastqInputs()
                input_bam_ch = handleBamInputs()
            } 
            else {
                error "Invalid file_types parameter. Must be one of: fastq, bam, bam_and_fastq"
            }
            
            //If there are any BAM inputs, convert them to FASTQ:
            converted_fastq = Channel.empty()
            if (input_bam_ch) {
                converted_fastq = Bam2Fastq(input_bam_ch)
            }

            // Combine BAM and FASTQ channels if both are provided
            def combined_fastq = null
            if (params.file_types == "fastq") {
                combined_fastq = input_fastq_ch
            }
            else if (params.file_types == "bam") {
                combined_fastq = converted_fastq
            }
            else if (params.file_types == "bam_and_fastq") {
                combined_fastq = input_fastq_ch.mix(converted_fastq)
            }           

            // Now, combined_fastq is a channel of tuples (name, technology, seqmode, fastq_file)
            // Duplicating the input channel for processes that consume it
            def (for_FASTQC, for_TrimAdaptors) = combined_fastq.multiMap { it ->
                    for_FASTQC: it
                    for_TrimAdaptors: it
                }

            // Run FASTQC on the combined input files
            fastqc_results = FASTQC(for_FASTQC)

            // Run TrimAdaptors on the combined inputs
            trimmed_fastq = TrimAdaptors(for_TrimAdaptors)
            // Run Bowtie2Alignment on trimmed FASTQ files
            alignment_outputs = Bowtie2Alignment(trimmed_fastq, bw2index.first())

            // Convert SAM to sorted BAM
            bam_outputs = Sam2Bam(alignment_outputs)

            // Removing duplicates using picard
            dedup_bam = rmdup(bam_outputs)

            // Sorting and indexing the dedup bam files
            index_rmdup(dedup_bam)

            // Collecting aligned BAMs and BAIs in a folder for ChromHMM
            Channel.value(file(params.outdir + '/Aligned_dedup_BAMS/'))
                   .set { aligned_dedup_indexed_BAMs_BAIs }

            // Channel for the cellmarkfiletable
            cellmarkfiletable = Channel.fromPath(params.cellmarkfiletable) 
            // Channel for the chrlengthfile
            chrlengthfile = Channel.fromPath(params.chrlengthfile)

            // ChromHMM BinarizeBam to binarize all files
            binarized_bams = binarize(aligned_dedup_indexed_BAMs_BAIs, cellmarkfiletable, chrlengthfile)

            // Parse the "min..max" string provided via the parameter (e.g. --numstates 2..50) 
            def (min, max) = params.numstates.split(/\.\./)*.toInteger()
            // Create a channel with the max value
            max_numstates = Channel.value(max)
            // Create a channel from the range between min and max (inclusive)
            numstates_ch = Channel.from(min..max)
            // Combining stuff needed for learnmodel in one channel 
            initial_channel = binarized_bams.combine(assembly_name)
            for_learnmodel = initial_channel.combine(numstates_ch)
            
            // ChromHMM learnmodel for ${numstates} number of states
            learned_model = learnmodel(for_learnmodel)

            // Running comparemodel for max states mentioned
           comparemodels(learnmodel.out.emissions.collect(), max_numstates)

            // Preparing channel for overlapenrichment
            feature_dir = Channel.fromPath(params.feature_dir)
            for_overlap = learnmodel.out.segments.map { file -> 
                                                        def name = file.name.replaceAll(/_segments\.bed\.gz$/, '')
                                                        tuple (name, file) }
                    .combine(feature_dir)

            // ChromHMM overlapenrichment on all segments files
            overlapenrichment(for_overlap)

            // Preparing channels for neighbourhood enrichment analysis
            TSS_ch = Channel.fromPath(params.TSS)
            TTS_ch = Channel.fromPath(params.TTS)

            // Combining segments files TSS and TTS bed files into a single tuple
            for_NE = learned_model[1].map { file ->
                                            def name = file.name.replaceAll(/_segments\.bed\.gz$/, '')
                                            tuple(name, file) }
                                        .combine(TSS_ch)
                                        .combine(TTS_ch)

            // Running Neighbourhood Enrichment analysis on all states models
            NeighbourhoodEnrichment(for_NE)

            // Read one or more files for marks/vars to be included or not
            def eval_subset_order = params.eval_subset_order instanceof String ? 
            params.eval_subset_order.split(',') : params.eval_subset_order
            Channel.fromPath(eval_subset_order.collect{ it.trim() })
                    .ifEmpty { error "No files found for input: ${params.input}" }
                    .map { file ->
                            def name = file.name.replaceAll(/\.txt$/, '') 
                            def content = file.text.trim()
                            tuple(name, content)}
                    .set { eval_subset_ch }

            segments_dir_ch = Channel.fromPath("$params.outdir/chromhmm/learnmodel/")

            // Joining model and segments channel by correct model number 
            learnmodel.out.allmodels.combine(segments_dir_ch)
                                    .combine(binarized_bams)
                                    .combine(eval_subset_ch)
                                    .set {eval_subset_input_ch}

            // Evaluating subset for user provided binary order/orders
           EvalSubset(eval_subset_input_ch) 
}

// Workflow completion handler
workflow.onComplete {
    log.info"""
    ========================================================================
    Workflow completed at: ${workflow.complete}
    Duration            : ${workflow.duration}
    Success            : ${workflow.success}
    workDir            : ${workflow.workDir}
    exit status        : ${workflow.exitStatus}
    ========================================================================
    """
}
