#!/usr/bin/env nextflow

// Enable DSL2 syntax (modern Nextflow with explicit workflow blocks)
nextflow.enable.dsl=2

// =========================================================================================
// IMPORT PROCESS MODULES
// =========================================================================================
// Modular architecture: Each process in separate file for maintainability

include { VALIDATE_INPUT }           	from './modules/validate_input.nf'
include { FASTQC as FASTQC_RAW }     	from './modules/fastqc.nf'
include { FASTQC as FASTQC_TRIMMED }	from './modules/fastqc.nf'     // Reusable process with alias
include { CELLRANGER_COUNT }			from './modules/cellranger_count.nf' 
include { MULTIQC }                  	from './modules/multiqc.nf'
//include { TEST_INDEX }           		from './modules/test_index.nf'  // Debugging utility

// =========================================================================================
// PIPELINE HEADER
// =========================================================================================
// Displays configuration for user verification before execution

log.info """
    ===========================================
    SCRNA-SEQ PIPELINE
    ===========================================
    Project          	: ${params.project}
    Species          	: ${params.species}
    Genome           	: ${params.genome_version}
    Fasta File       	: ${params.ref_fasta()}
    GTF File         	: ${params.ref_gtf()}
    BED File         	: ${params.ref_bed()}
    Housekeeping     	: ${params.housekeeping_bed()}
    
    PATHS:
    Reference Dir    	: ${params.ref_dir()}
    STAR Index       	: ${params.star_index_dir()}
    Salmon Index     	: ${params.salmon_index_dir()}
    Project Dir      	: ${params.proj_dir()}
    Input (FastQ)    	: ${params.fastq_dir()}
    Input (Raw FastQ)	: ${params.raw_fastq_dir()}
    Output (FastQC)  	: ${params.fastqc_dir()}   
    Output (CellRanger)	: ${params.cellranger_dir()}
    Output (MultiQC) 	: ${params.multiqc_dir()}
    Logs             	: ${params.log_dir()}
    ===========================================
    """

// =========================================================================================
// MAIN WORKFLOW
// =========================================================================================
// Data flows through processes via channels
// Nextflow automatically parallelizes based on channel cardinality

workflow {
    
    // =====================================================================================
    // STEP 1: VALIDATE INPUT FASTQ FILES
    // =====================================================================================
    // Validates naming conventions and creates organized sample channels
    
    VALIDATE_INPUT(params.raw_fastq_dir())
    mode_ch         = VALIDATE_INPUT.out.mode.collect()     // [mode] → collected for RSeQC
    sample_fastq_ch = VALIDATE_INPUT.out.samples            // [sample_id, [R1, R2]]
    
    // =====================================================================================
    // STEP 2: QUALITY CONTROL ON RAW READS
    // =====================================================================================
    // FastQC analyzes read quality before processing
    
    // Add read_type tag to channel: [sample_id, [fastqs], "raw"]
    fastqc_ch = sample_fastq_ch.map { sample_id, fastqs -> tuple(sample_id, fastqs, "raw") }
    FASTQC_RAW(fastqc_ch)
    
    // =====================================================================================
    // STEP 3: CELLRANGER COUNT
    // =====================================================================================   
    // CRITICAL: Pre-join arguments to prevent cache invalidation
    // If params.CELLRANGER_ARGS().join(' ') called inside process → hash changes → resume fails
    cellranger_args = params.CELLRANGER_ARGS().join(' ')
    CELLRANGER_COUNT(sample_ch, cellranger_args)    

    // =====================================================================================
    // STEP 4: AGGREGATE QC REPORTS (MULTIQC)
    // =====================================================================================
    // Combines all QC outputs into single interactive HTML report
    
    multiqc_ch = Channel.empty()
        .mix(FASTQC_RAW.out.fastqc_zip)                     // FastQC reports
        .mix(SALMON_QUANT.out.salmon_quant_dir.map { it[1] })   // Salmon dirs (extract from tuple)
        .mix(STAR_ALIGN.out.gene_counts)                    // ReadsPerGene.out.tab
        .mix(STAR_ALIGN.out.sj_tab)                         // Splice junctions
        .mix(STAR_ALIGN.out.star_log)                       // Alignment stats
        .mix(RSEQC.out.rseqc_logs)                          // RSeQC outputs
        .collect()                                          // Wait for all samples
    
    // CRITICAL: Convert closures to string and pass into process to  prevent cache invalidation
    multiqc_title = params.multiqc_titlename()
    multiqc_file  = params.multiqc_filename()	
    MULTIQC(multiqc_ch, multiqc_title, multiqc_file)
}

// =========================================================================================
// PIPELINE OVERVIEW
// =========================================================================================
//
// This pipeline performs comprehensive RNA-seq analysis:
//   1. Input validation and quality control (FastQC)
//   2. Reference genome indexing (STAR, Salmon, BED conversion)
//   3. Transcript quantification (Salmon - fast, alignment-free)
//   4. Read alignment (STAR - splice-aware, generates BAM)
//   5. Alignment QC (RSeQC - detects biases and issues)
//   6. Report aggregation (MultiQC - single HTML report)
//
// OUTPUT STRUCTURE:
// ${proj_dir}/
// ├── 01.FastQ/raw/                      # Raw input FASTQs
// ├── 02.FastQC/raw/                     # QC on raw reads  
// ├── 03.Salmon/                         # Transcript quantification
// │   ├── Sample1/                       # Full Salmon output
// │   └── quant_files/                   # Collected quant.sf files
// ├── 04.STAR/                           # Alignment outputs
// │   ├── gene_counts/                   # ReadsPerGene.out.tab files
// │   ├── splice_junction/               # SJ.out.tab files
// │   ├── alignment_stats/               # Log.final.out files
// │   └── bam			                  # BAM + BAI files
// ├── 05.RSEQC/                          # Organized by analysis type
// │   ├── 01_read_distribution/
// │   ├── 02_inner_distance/
// │   ├── 03_junction_annotation/
// │   └── ... (09 subdirectories total)
// ├── 06.MultiQC/
// │   ├── Project_MultiQC_Report.html
// │   └── Project_MultiQC_Report_data/
// └── 07.Logs/                           # All error logs + reports
//     ├── trace.txt
//     ├── report.html
//     └── timeline.html
//
// DOWNSTREAM ANALYSIS:
//   Differential Expression:
//     - Option 1: Use STAR gene counts with DESeq2/edgeR directly
//     - Option 2: Import Salmon quant.sf files with tximport → DESeq2/edgeR
//   
//   Visualization:
//     - Load BAM files into IGV for gene-level inspection
//     - Use RSeQC or deepTools for coverage plots
//   
//   Splice Analysis:
//     - Use STAR SJ.out.tab files with rMATS or LeafCutter
//
// QUALITY CONTROL THRESHOLDS (check MultiQC report):
//   - Mapping rate: >70% (good), 60-70% (acceptable), <60% (investigate)
//   - Duplication: <50% (good), 50-70% (acceptable), >70% (low complexity)
//   - Read distribution: >50% CDS exons, <10% introns, <5% intergenic
//   - 3' bias: <3 (good), 3-5 (acceptable), >5 (degraded RNA)
//
// RESUMING FAILED RUNS:
//   Nextflow caches completed processes in work/ directory
//   To resume: bash run_nextflow.sh (script includes -resume flag)
//   Cache is invalidated if: process code changes, input files change, or work/ deleted
//
// CLEANING UP:
//   After successful run and verification:
//     rm -rf ${work_dir}/*              # Free disk space
//   Warning: Cannot resume after deleting work directory!
//
// COMMON ISSUES:
//   "No such file or directory":
//     → Check paths in project_info.yaml
//     → Verify Singularity bind mounts (run_nextflow.sh prints mappings)
//   
//   Out of memory:
//     → Increase memory in process labels (nextflow.config)
//     → Reduce parallel jobs with maxForks in nextflow.config
//   
//   Process won't resume:
//     → Check if params.FUNCTION().join(' ') called inside process (should be in workflow)
//     → Verify work/ directory not deleted
//     → Review .nextflow.log for cache invalidation reasons
//
// EXTENDING THE PIPELINE:
//   1. Create new process module in modules/ directory
//   2. Add include statement above (top of file)
//   3. Call process in workflow block (add to appropriate step)
//   4. Configure publishDir in nextflow.config
//   5. Add outputs to multiqc_ch if process generates QC data
//
// For detailed documentation on each process:
//   docs/validate_input.md, docs/star_align.md, docs/salmon_quant.md, etc.
//
// =========================================================================================
