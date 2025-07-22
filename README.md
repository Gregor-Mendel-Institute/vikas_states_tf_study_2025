Github repository for Chromatin states pridiction nextflow pipeline, preprocessed data, and R code to reproduce manuscript figures.

# 1. chromatin_states_discovery_nextflow
DSL2 nextflow pipeline for chromatin state pridiction using chromHMM. The pipeline is designed to take raw ChIP-seq data as input and performs the following tasks: 
1. Pull genome index from singularity or built custom genome index using bowtie2 build
2. Convert BAM files to FASTQ if provided
3. Perform initial QC using FASTQC
4. Trims reads using trim_galore
5. Perform final QC using FASTQC
6. Alignment to reference genome (TAIR10, MpTak1v5.1r2 using singularity container or any custom reference genome) using bowtie2
7. Removal of duplicates using picard
8. Chromatin state pridiction using binarizeBam and learnModel utilities of chromHMM
9. Evaluate enrichment of each state in genomic features using overlapenrichment utility of chromHMM
10. Perform Neighbourhood Enrichment analysis with specific anchor positions using NeighbourhoodEnrichment utility of chromHMM
11. Perform model pridiction using sub-epigenomes using EvalSubset utility of chromHMM
12. Generate BigWig files using the bamCompare tool of deeptools
13. ComputeMatrix and generate heatmap and metaplots of ChIP signal on specific regions using ComputeMatrix and plotHeatmap tools of deeptools

## Raw Input Files
The pipeline can accept:
1. Raw BAM files of SR, PE or mixed sequencing mode. (Where mixed means some files of SR and some of PE sequencing mode)
2. Raw FASTQ files of SR, PE or mixed sequencing mode.
3. Raw BAM and raw FASTQ files togather with any combination of sequencing modes.
The pipeline can also recognize the type of experiment (ChIP-seq or CUT&RUN) and process it accordingly.

To input raw BAM files use the glob patterns like this:
```
--bams "*.bam"
```
If one wishes to input raw FASTQ files, input the raw files using glob pattern like this:
For SR FASTQ files:
```
--fqs_sr "*.fastq"
```
For PE FASTQ files:
```
--fqs_pe "*_{1,2}.fastq"
```
Do not forget to put your glob pattern in **""** as mentioned above else nextflow will only read a single file instead of all files matching the pattern! 

Provide the information of the input file types to be used in the pipeline using the parameter **file_types**. There are three possible input file_types:
1. In case of only raw FASTQ files
```
--file_types fastq
```
2. In case of only raw BAM files
```
--file_types bam
```
3. In case of both BAM and FASTQ files
```
--file_types bam_and_fastq
```

In case of file_types **fastq** or **bam_and_fastq** an additional parameter of **fastq_file_type** need to be provided which again has three possibilities:
1. In case of only SR FASTQ files
```
--fastq_file_type Only_SR
```
2. In case of only PE FASTQ files
```
--fastq_file_type Only_PE
```
3. In case of both SR and PE files
```
--fastq_file_type SR_and_PE
```

## Input sample information file
The information about the input files should be provided in a tab seperated file with three columns. First column should contain the file name **(without the file extension)** , the second column should contain the name of sequencing technology ('ChIP' for ChIP-seq files and 'CNR' for cut&run files) and the third column should contain the sequencing mode (For paired end files use: **PE** and for single end file use: **SR**) 
One eg. of file layout is shown here:

file1            ChIP          SR

file2            ChIP          SR

file3            CNR           PE

file4            CNR           PE

(Note that there are no file extensions in this file even though the file is originally called file1.bam or file1.fastq)

This file can now be parsed to the pipeline using the parameter 
``` 
--bam_samples bam_samples.tab
```

If you wish to parse fastq files then create a similar file with sample information and parse it like this
```
--fq_samples fq_samples.tab
```

## CELLMARKFILETABLE for ChromHMM
