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
ChromHMM requires a tab seperated cellmarkfiletable which is a file containing atleast 3 columns with cell type, mark and file name. An additional fourth column for the control file can also be provided for normalization (See ChromHMM manual). The file layout should be like this:

MP      H2A.X  test.dedup.sorted.bam control.dedup.sorted.bam

MP      H2A.W  test.dedup.sorted.bam control.dedup.sorted.bam

MP      H2A.Z   test.dedup.sorted.bam control.dedup.sorted.bam

Parse this file to the pipeline with following parameter
``` 
--cellmarkfiletable cellmarkfiletable.tab
```

## CHROMOSOME LENGTH FILE

A chromosome length file with two columns where the first column contains the chromosome name and second column with it's length with following layout

chr1    30584173

chr2    29643427

chr3    27142341

Make sure that the chromosome names in this file are same as the assembly that it's going to be aligned to. For Arabidopsis, remove the 'chr' from this file!
Parse this file to the pipeline with following parameter
``` 
--chrlengthfile chrlengthfile.txt 
```

## EVAL-SUBSET BINARY STRING FILES
To evaluate the sub-epigenome models using only a subset of marks using EvalSubset utility of chromHMM. Provide a comma seperated list of text files each containing a binary string where 0 is to ignore the mark and 1 is to keep it for the model pridiction. See chromHMM manual for more information. Example usage:
```
--eval_subset_order eval/all_with_all.txt,eval/only_marks.txt,eval/only_vars.txt
```
Where files look like this:
```
eval/all_with_all.txt --> 111111111111111
eval/only_marks.txt --> 000011111111111
eval/only_vars.txt --> 111100000000000
```

## Parameters
### NEXTFLOW PARAMETERS
- **-profile:** slurm or local.
- **-w:** Working directory where all the temperory files will be store. It should **ALWAYS** be in the /scratch-cbe.
- **--outdir:** the name of the output folder. [ default: results ]
- **-resume:** To resume the pipeline run interrupted for some reason.    
- **-bg:** To run the pipeline in the background.
  
### REFERENCE GENOME PARAMETERS 
- **--reference_fasta:** FASTA file for custom reference genome.
- **--reference_genome:** tair10, mptak(MpTak1v5.1), or grcm38(mm10) for singularity containers. Either choose an existing singularity container else build custom genome index by providing FASTA file for refernce genome using --reference_fasta parameter.
- **--index_basename:** Variable foor the genome assembly name
  
### INPUT FILE PARAMETERS
- **--file_types:**  bam or bam_and_fastq or fastq.
- **--bam_samples:** If --file_types:  bam or bam_and_fastq then proveide a tab seperated file with three columns (1. File name without extension 2. Sequencing mode (SR or PE) 3. Sequencing technology (ChIP or CNR) of all bam files.
- **--fq_samples:** If --file_types:  fastq or bam_and_fastq then proveide a tab seperated file with three columns (1. File name without extension 2. Sequencing technology (ChIP or CNR) and 3. Sequencing mode (SR or PE) of all fastq files.
- **--fastq_file_type:** If --file_types: bam_and_fastq or fastq then this parameter can take three values: Only_SR, Only_PE or SR_and_PE.
- **--bams:** If --file_types:  bam or bam_and_fastq then with this parameter provide all the bam files using wild card eg. "*.bam". Don't forget to add "" otherwise it will only read one file.
- **--fqs_sr:** If --file_types: bam_and_fastq or fastq and --fastq_file_type Only_SR then with this parameter provide all the SR fastq files using wild card eg. "*.fastq". Don't forget to add "" otherwise it will only read one file.
- **--fqs_pe:** --file_types: bam_and_fastq or fastq and --fastq_file_type Only_PE or SR_and_PE then with this parameter provide all the PE fastq files using wild card eg. "*_{1,2}.fastq". Don't forget to add "" otherwise it will only read one file.
  
### chromHMM PARAMETERS
- **--cellmarkfiletable:** Tab seperated file as mentioned above with three or four columns. See ChromHMM manual for more information.
- **--chrlengthfile:** Tab seperated file as mentioned above with two columns.
- **--numstates:** Provide a range of number of states to be built by the learnmodel. Parse the "min..max" string provided via the parameter (e.g. --numstates 2..50)
- **--feature_dir:** A directory with BED files for genomic features like gene.bed, five_prime_utr.bed etc.
- **--TSS:** A BED file containing the position of TSS used in neighbourhood-enrichment calculation by chromHMM. 
- **--TTS:** A BED file containing the position of TTS used in neighbourhood-enrichment calculation by chromHMM.
- **--eval_subset_order:** Provide comma seperated list of text files containing a binary string each, referencing to which marks to keep (1) and which marks to ignore (0) for the sub-epigenomes. See chromHMM manual for more information.

## OUTPUT DIRECTORIES AND FILES
The pipeline should create the following main directories within the mentioned output directory provided by --outdir parameter. 

1. **Genome Index Directory**: If the custom genome index is build the genome is saved in a directory named as the --index_basename parameter. 
2. **First_Fastqc_Reports**: This directory shall contain the Fastqc reports of all the input files before any read trimming.
3. **FilterFastq**: This directory shall contain the fastqc reports after adaptor trimming.
4. **Alignment_logs**: This directory shall contain the alignment logs of each file aligned to reference genome using bowtie2.
5. **Aligned_dedup_BAMS**: This directory shall contain the duplicate removal logs from picard run on each file.
9. **chromhmm**: This directory contains sub-directories containing chromHMM output files from different utilities. The sub-directories that shall be created within this folder should be the following:
    **A. chromhmm/binarized_files**: This directory shall contain the binarized data files created by chromHMM for model learning.
    **B. chromhmm/learnmodel**: This directory shall contain the chromatin state models emitted by chromhmm learnModel.
    **C. chromhmm/comparemodels**: This directory shall contain the output of chromHMM comparemodels.
    **D. chromhmm/overlapenrichment**: This directory contains the enrichment of each state in the genomic features calculated by overlapenrichment utility of chromHMM.
    **E. chromhmm/Neighbourhood_Enrichment**: This directory contains the enrichment of each state with respect to TSS and TTS calculated by neighbourhood utility of chromHMM. The results of the two anchor positions TSS and TTS are placed in two sub-directories within this directory named **chromhmm/Neighbourhood_Enrichment/TSS** and **chromhmm/Neighbourhood_Enrichment/TTS**.
    **F. chromhmm/evalSubset**: This directory contains as many sub-directories with the same name as many are there files provided for the parameter --eval_subset_order each containing the sub-epigenome models.

## RECOMMENDED SETUP

#### -1. Add Hinkskalle to singularity (one time only)
#### 0. Start a new tmux session or attach to an existing (optional but useful)

See e.g. https://tmuxcheatsheet.com/ 

#### 1. Clone the repository

From  e.g. *your user folder in the lab folder* do:

```
git clone https://github.com/Gregor-Mendel-Institute/vikas_states_tf_study_2025.git
```
you might be asked for user and password, use the ones from your **forskalle** account.

 #### 2. Create a folder with links to the bam/fastq

e.g:
```
mkdir -p bams
ln -s /groups/berger/Raw/demultiplexed/97009_AATGATAGGACCACCT_CDR8FANXX_1_20190821B_20190821.bam bams/
ln -s /groups/berger/Raw/demultiplexed/97010_TTGGCGCGTGGCTAGG_CDR8FANXX_1_20190821B_20190821.fastq fastqs/
ln -s /where/ever/you/have/your/file_1.fastq PE_fastqs/
ln -s /where/ever/you/have/your/file_2.fastq PE_fastqs/
```    
As long as you use links and do not actually copy the files, this folder can be created in your user folder on the lab folder, e.g. in the Chromatin_state_discovery_nf folder that you cloned or in the folder above. 
**If you have fastq files of sequencing type PE and SE both, then it's a good idea to keep then in seperate folder sothat you can input SR files with wild card * (eg. bams/*.bam) and PE file pairs with wild card *_{1,2}.fastq. Keeping them in the same folder might create a clash so be careful here!**

#### 3. Then load the following modules:

```
module load nextflow/21.04.1  #or higher if you prefer
singularity remote use hinkskalle
```  
**NOTE** if you are in a tmux session were this was already done, this step can be skipped!

### 4.  Go to the vikas_states_tf_study_2025 folder.

```
cd vikas_states_tf_study_2025/
```
## How to run 
 
 **Make sure you have:**
 1. completed the setup steps 
 2. that you are on the cluster with the correct modules loaded
 3. moved into the vikas_states_tf_study_2025/ folder

To use the pipeline with the only BAM files  
```
nextflow run main.nf -with-dag flowchart.html -profile slurm --file_types "bam" --bam_samples "sample_Info.tab" --bams "/groups/berger/user/vikas.shukla/MP_states/bams/*.bam" --reference_fasta "/groups/berger/user/vikas.shukla/v7.1_MP_states/Genome/Fasta/MpTak1_v7.1.fa.gz" --cellmarkfiletable "/groups/berger/user/vikas.shukla/v7.1_MP_states/Chromatin_state_discovery_nf/cellmarkfiletable.tab" --chrlengthfile "/groups/berger/user/vikas.shukla/v7.1_MP_states/Chromatin_state_discovery_nf/MPTak1_v7_1_chrlengths.txt" --index_basename MpTak1_v7.1 --numstates 2..50 --feature_dir "/groups/berger/user/vikas.shukla/v7.1_MP_states/Genome/Annotations/Features/" --TSS /groups/berger/user/vikas.shukla/v7.1_MP_states/Genome/Annotations/TSS.bed --TTS /groups/berger/user/vikas.shukla/v7.1_MP_states/Genome/Annotations/TES.bed --eval_subset_order eval/all_with_all.txt,eval/only_marks.txt,eval/only_vars.txt --outdir "MpTak1v7.1_states_06052025/" -w "/scratch-cbe/users/vikas.shukla/MP_states_v7.1_06052025/" 
```

The -w sets the path to where the pipeline should be run and where the intermediate files will be stored. **It should be on the scratch and not in the lab folder!**    

To use the pipeline with the mixed file types (that is **BAM and FASTQ files*)
```
nextflow run main.nf -profile slurm --file_types bam_and_fastq --bam_samples "../bams/sample_Info.tab" --bams "../bams/*.bam" --fq_samples "fastq_samples.tab" --fastq_file_type "Only_PE" --fqs_pe "../bams/ATACseq/*.fastq" --reference_genome "mptak" --chrlengthfile mp_chrlengthfile.txt --cellmarkfiletable "mp_cellmarkfiletable.tab" --numstates 2..50 --featureDir "../Mp_genome/newFeatures/" --TSS "../Mp_genome/TSS.bed" --TTS "../Mp_genome/TTS.bed" --eval_subset_order eval/all_with_all.txt,eval/only_marks.txt,eval/only_vars.txt -w "/scratch-cbe/users/vikas.shukla/new_results/" --outdir "New_results/" 
```


## <a name="hinkskalle"></a>How to adding Hinkskalle to singularity

As the pipeline uses containers from the ngs registry (aka "Hinkskalle"), you need to do the following (on the CBE cluster):

```
singularity remote add --no-login hinkskalle singularity.ngs.vbcf.ac.at
singularity remote use hinkskalle
```

**It's enough if you do this once.** The settings will then be saved a folder (.singularity) in your home directory and you do not need to worry about it any more. **UNLESS** you later change the remote registry ( singularity remote use <another registry> ). Then you have to re-run the last line before you start the pipeline.

**NOTE:** If you are using another cluster (e.g. outside of VBC ), you need to make sure Singularity (version 3.4 or higher) is installed and loaded BEFORE you add the hinklist registry. 


##  <a name="tmux"></a>Short intro to tmux

tmux is a very convent tool to use on the cluster. It has many benefits but some of the most obvious are:

1. It lets you set up an environment on the cluster that you can move in an out from as much you want. E.g. if you have a tmux session called ChIPseq_NF you can load the modules needed to run the pipeline (nextflow and singularity) **once** and then by re-attaching to this session the modules are already loaded.
2. Scripts that are running in the session will continue to run even if you logout from the cluster. 

See also e.g.  https://en.wikipedia.org/wiki/Tmux/

To start a new session type:

    tmux new -s vikas_states_tf_study_2025

This will create a tmux session with the name vikas_states_tf_study_2025 and attach you to it. 
