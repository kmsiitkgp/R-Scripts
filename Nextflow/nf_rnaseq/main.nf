#!/usr/bin/env nextflow

// =========================================================================================
// RNA-SEQ ANALYSIS PIPELINE
// =========================================================================================
// Purpose: Comprehensive RNA-seq analysis pipeline for differential gene expression
//
// Pipeline steps:
//   1. Input validation and quality control
//   2. Reference genome indexing (STAR, SALMON)
//   3. Read alignment (STAR)
//   4. Transcript quantification (SALMON)
//   5. Alignment quality control (RSeQC)
//   6. Report aggregation (MultiQC)
//
// Supported modes:
//   - Single-end (SE) sequencing
//   - Paired-end (PE) sequencing
//   - Tumor/Normal sample pairs
//
// Requirements:
//   - Nextflow >= 21.0
//   - Singularity containers (configured in nextflow.config)
//   - Reference genome files (FASTA + GTF)
//   - FASTQ files in expected directory structure
//
// Usage:
//	 module load nextflow/24.10.5
//	 module load singularity-apptainer/1.1.8
//   nextflow run main.nf -profile <profile_name> [options]
//
// Common options:
//   -profile <name>    : Select project profile (required)
//   -resume            : Resume from last successful step
//   -with-dag dag.png  : Generate workflow diagram
//   -with-report       : Generate execution report
//   -with-timeline     : Generate timeline visualization
//
// Example:
//   nextflow run main.nf -profile xinyi -resume
//
// =========================================================================================

// Enable DSL2 syntax (modern Nextflow syntax with explicit workflow blocks)
nextflow.enable.dsl=2

// =========================================================================================
// IMPORT PROCESS MODULES
// =========================================================================================
// Each module contains a single process definition
// This modular approach improves code organization and reusability
// All processes are defined in separate files under the ./modules/ directory

include { VALIDATE_INPUT }           from './modules/validate.nf'       // Validates FASTQ naming and structure
include { FASTQC as FASTQC_RAW }     from './modules/fastqc.nf'         // QC on raw reads
include { FASTQC as FASTQC_TRIMMED } from './modules/fastqc.nf'         // QC on trimmed reads (reuses same process)
include { PREP_REFERENCE }           from './modules/prep_reference.nf' // Index genome for STAR and SALMON
include { SALMON_QUANT }             from './modules/salmon_quant.nf'   // Transcript quantification
include { STAR_ALIGN }               from './modules/star_align.nf'     // Genome alignment
include { RSEQC }                    from './modules/rseqc.nf'          // Alignment QC
include { MULTIQC }                  from './modules/multiqc.nf'        // Aggregate reports

// =========================================================================================
// PROFILE VALIDATION
// =========================================================================================
// Ensures user has selected a valid project profile before running
// Prevents accidental runs with default "test" parameters
// This check runs before any processing begins, saving time and resources

if ( params.species == 'test' || params.project == 'test' ) {
    error """
❌ No project profile selected!

You must run this pipeline with a project profile, for example:

    nextflow run main.nf -profile xinyi

Available profiles:
  Human projects:
  - manish_arc     : Manish_22RV1_ARCaPM project
  - manish_xeno    : Manish_22RV1_Xenograft project
  - manish_xeno2   : Manish_22RV1_Xenograft2 project
  - vera           : Vera project
  - xinyi          : Xinyi project

  Mouse projects:
  - sandrine       : Sandrine_Quad project
  - supriya        : Supriya_Quad project

To add a new profile, edit the profiles section in nextflow.config
"""
}

// =========================================================================================
// PIPELINE HEADER
// =========================================================================================
// Displays key parameters for user verification
// Helps ensure correct genome/project is being used before processing begins
// This information is printed to the console when the pipeline starts

log.info """
    ===========================================
    RNA-SEQ PIPELINE
    ===========================================
    Project : ${params.project}
    Species : ${params.species}
    Genome  : ${params.genome_version}
    ===========================================
    """

// =========================================================================================
// MAIN WORKFLOW
// =========================================================================================
// Defines the execution flow and data dependencies between processes
// Nextflow automatically parallelizes based on channel cardinality
// Each process runs as soon as its input dependencies are satisfied
//
// Workflow execution order:
//   1. VALIDATE_INPUT (runs first, requires no dependencies)
//   2. FASTQC_RAW and PREP_REFERENCE (run in parallel after step 1)
//   3. SALMON_QUANT and STAR_ALIGN (run in parallel after step 2)
//   4. RSEQC (runs after STAR_ALIGN completes)
//   5. MULTIQC (runs last, after all QC processes complete)
	
workflow {

	// =====================================================================================
	// STEP 1: VALIDATE INPUT FASTQ FILES
	// =====================================================================================
	// Purpose: Ensure FASTQ files follow naming conventions and detect sequencing mode
	// 
	// Input:
	//   - Directory path (string) from params.raw_fastq_dir()
	//   - Expected FASTQ naming: SampleID_R1.fastq.gz, SampleID_R2.fastq.gz
	//
	// Output:
	//   - mode: "SINGLE_END" or "PAIRED_END" (based on detected files)
	//   - samples: Channel of [sample_id, [fastq_files]] tuples
	//   - Additional metadata: sample tags, counts, etc.
	//
	// Validation checks:
	//   - Files exist and are readable
	//   - Naming conventions are followed
	//   - Consistent sequencing mode across all samples
	//   - Paired files are properly matched
	//
	// Note: Workflow can accept strings, but processes require channels
	//       VALIDATE_INPUT creates channels from the input directory
	// =====================================================================================
	
    VALIDATE_INPUT(params.raw_fastq_dir())
	
	// Extract outputs from VALIDATE_INPUT for use in downstream processes
	mode            = VALIDATE_INPUT.out.mode.collect()  // "SINGLE_END" or "PAIRED_END"
	sample_fastq_ch = VALIDATE_INPUT.out.samples         // Channel: [sample_id, [fastq_files]]
	
	// =====================================================================================
	// STEP 2: QUALITY CONTROL ON RAW READS
	// =====================================================================================
	// Purpose: Assess raw read quality before any processing
	// Runs FastQC on all input FASTQ files to identify potential issues
	//
	// Quality metrics assessed:
	//   - Per-base sequence quality
	//   - Per-sequence quality scores
	//   - Per-base sequence content
	//   - GC content distribution
	//   - Sequence length distribution
	//   - Sequence duplication levels
	//   - Overrepresented sequences
	//   - Adapter content
	//
	// Data transformation:
	//   Input:  [sample_id, [R1.fq.gz]] or [sample_id, [R1.fq.gz, R2.fq.gz]]
	//   Output: [sample_id, [R1.fq.gz], "raw"] or [sample_id, [R1.fq.gz, R2.fq.gz], "raw"]
	//
	// The "raw" tag is used by FASTQC for output directory organization
	// This allows separating raw vs. trimmed read QC reports
	// Multiple FASTQC processes run in parallel (one per sample)
	// =====================================================================================
	
	fastqc_ch = sample_fastq_ch.map { sample_id, fastqs -> tuple(sample_id, fastqs, "raw") }
	FASTQC_RAW(fastqc_ch)
	
	// =====================================================================================
	// STEP 3: PREPARE REFERENCE GENOME INDEXES
	// =====================================================================================
	// Purpose: Build genome indexes for alignment and quantification
	// This is a computationally expensive step but only needs to run once per genome
	//
	// Creates three main outputs:
	//   - STAR index: For splice-aware genome alignment
	//     - Contains SA (suffix array), SAindex, Genome files
	//     - Enables fast seed matching and extension
	//   - SALMON index: For fast transcript quantification
	//     - Contains k-mer based index for pseudo-alignment
	//     - Much smaller than genome index
	//   - BED file: For RSeQC quality control
	//     - Gene annotation in BED12 format
	//     - Used for read distribution analysis
	//
	// Indexing parameters:
	//   - STAR uses parameters from params.STAR_ARGS
	//   - SALMON uses k-mer size optimized for transcript quantification
	//
	// Note: Indexes are published to ref_dir (not project directory)
	//       This allows reuse across multiple pipeline runs and projects
	//       Once created, indexes can be reused indefinitely unless genome changes
	// =====================================================================================
	
	// Create value channels for reference files
	// Channel.value() creates a singleton channel (only one item, reused by all processes)
	// checkIfExists: true causes pipeline to fail immediately if files are missing
	// This prevents wasting time on processing before discovering missing references
	ref_fasta_ch = Channel.value(file(params.ref_fasta[params.species], checkIfExists: true))
    ref_gtf_ch   = Channel.value(file(params.ref_gtf[params.species], checkIfExists: true))
	
	// Run reference preparation process
	PREP_REFERENCE(ref_fasta_ch, ref_gtf_ch)
	
	// Extract output channels for use in downstream processes
	// Need to use .collect() on single item channels because in the next step
	// we combine multi-item channels (sample_fastq_ch) with single-item channels
	// Nextflow will throw an error if the single-item channel was not built using .collect()
	// .collect() ensures proper channel behavior when mixing cardinalities
	star_index_ch   = PREP_REFERENCE.out.star_index_dir.collect()   // STAR genome index directory
	salmon_index_ch = PREP_REFERENCE.out.salmon_index_dir.collect() // SALMON transcript index directory
	ref_bed_ch      = PREP_REFERENCE.out.ref_bed.collect()          // BED annotation file
	
	// =====================================================================================
	// STEP 4: QUANTIFY TRANSCRIPTS USING SALMON
	// =====================================================================================
	// Purpose: Fast, alignment-free transcript abundance estimation
	// Uses quasi-mapping and streaming fragment assignment for speed
	//
	// Algorithm overview:
	//   1. Map reads to transcripts using k-mer based pseudo-alignment
	//   2. Assign multi-mapping reads using expectation-maximization (EM)
	//   3. Estimate transcript abundances in TPM and read counts
	//   4. Apply bias corrections (GC, sequence, positional)
	//
	// Inputs:
	//   - sample_fastq_ch: [sample_id, [fastq_files]]
	//   - salmon_index_ch: Path to SALMON index directory
	//
	// Outputs:
	//   - quant.sf: Main quantification file (TPM, NumReads per transcript)
	//   - Sample-specific directories with auxiliary information
	//   - Mapping statistics and parameters used
	//
	// Advantages over alignment-based methods:
	//   - Much faster than alignment-based methods (~10-30x speedup)
	//   - Better handling of multi-mapping reads (uses EM algorithm)
	//   - Corrects for GC, sequence, and positional biases
	//   - Directly estimates transcript-level abundances
	//   - Lower memory requirements
	//
	// Note: Can run in parallel with STAR_ALIGN (independent processes)
	//       Both use the same input FASTQ files but serve different purposes
	// =====================================================================================
	
	SALMON_QUANT(sample_fastq_ch, salmon_index_ch)
	
	// =====================================================================================
	// STEP 5: ALIGN READS USING STAR
	// =====================================================================================
	// Purpose: Splice-aware alignment to reference genome
	// STAR (Spliced Transcripts Alignment to a Reference) is optimized for RNA-seq
	//
	// Algorithm overview:
	//   - Two-pass mode enabled for novel junction discovery
	//   - First pass: Map reads and discover splice junctions
	//   - Second pass: Re-map reads using discovered junctions
	//   - This improves mapping around novel splice sites
	//
	// Inputs:
	//   - sample_fastq_ch: [sample_id, [fastq_files]]
	//   - star_index_ch: Path to STAR index directory
	//
	// Outputs:
	//   - Coordinate-sorted BAM files with indexes (.bam, .bai)
	//     - Sorted by genomic coordinates for efficient access
	//     - Indexed for fast random access (required for IGV, RSeQC)
	//   - Gene count matrices (ReadsPerGene.out.tab)
	//     - Column 1: Gene ID
	//     - Column 2: Unstranded counts
	//     - Column 3: Counts for strand 1
	//     - Column 4: Counts for strand 2
	//   - Splice junction tables (SJ.out.tab)
	//     - Contains all detected junctions with supporting reads
	//   - Alignment statistics (Log.final.out)
	//     - Mapping rates, multi-mapping statistics, etc.
	//
	// Use cases:
	//   - Differential gene expression (use gene counts)
	//   - Visualization in IGV (use BAM files)
	//   - Quality control with RSeQC (use BAM files)
	//   - Splice variant analysis (use junction tables)
	//   - Novel junction discovery (two-pass mode)
	// =====================================================================================
	
	STAR_ALIGN(sample_fastq_ch, star_index_ch)
	
	// Extract BAM files with indexes for downstream QC
	// Channel contains: [sample_id, bam_file, bai_file]
	sample_bam_ch = STAR_ALIGN.out.bam_indexed
	
	// =====================================================================================
	// STEP 6: QUALITY CONTROL ON ALIGNMENTS
	// =====================================================================================
	// Purpose: Comprehensive QC analysis on aligned BAM files
	// RSeQC provides RNA-seq specific quality control metrics
	//
	// Inputs:
	//   - sample_bam_ch: [sample_id, bam, bai]
	//   - ref_bed_ch: Gene annotation in BED format
	//   - mode: "SINGLE_END" or "PAIRED_END"
	//
	// Analyses performed:
	//   - Read distribution across genomic features
	//     - Percentage of reads in CDS, UTRs, introns, intergenic regions
	//     - Helps identify library preparation issues
	//   - Gene body coverage uniformity (3' bias detection)
	//     - Plots coverage from 5' to 3' end of genes
	//     - Detects RNA degradation or library prep bias
	//   - Splice junction annotation
	//     - Known vs. novel junctions
	//     - Junction read support
	//   - Insert size distribution (PE only)
	//     - Fragment size distribution for paired-end data
	//     - Helps verify library preparation
	//   - Sequencing artifact profiling
	//     - Duplication rates
	//     - Read complexity
	//
	// Outputs:
	//   - PDF plots for visualization
	//   - Text/log files for MultiQC aggregation
	//   - R script files for custom plotting
	//
	// Quality thresholds to check:
	//   - >80% of reads should map to exons
	//   - 3' bias ratio should be <3
	//   - <10% of reads in intergenic regions
	// =====================================================================================
	
	RSEQC(sample_bam_ch, ref_bed_ch, mode)
	
	// =====================================================================================
	// STEP 7: AGGREGATE ALL QC REPORTS USING MULTIQC
	// =====================================================================================
	// Purpose: Combine QC reports from all tools into a single interactive HTML report
	// MultiQC automatically detects and parses output files from supported tools
	//
	// Strategy: Mix all QC outputs into a single channel, then collect them
	//
	// Data flow:
	//   1. Create empty channel as starting point
	//   2. Mix in outputs from each QC tool (FastQC, SALMON, STAR, RSeQC)
	//   3. Collect all items into a single list
	//   4. Pass complete list to MULTIQC
	//
	// Why use .collect()?
	//   - Waits for ALL samples to complete before running MultiQC
	//   - Ensures the report contains data from every sample
	//   - Creates a single MultiQC report instead of one per sample
	//   - Allows cross-sample comparisons and batch effect detection
	//
	// MultiQC features:
	//   - Interactive plots with sample filtering
	//   - Summary statistics across all samples
	//   - Outlier detection
	//   - Customizable report sections
	//   - Export to CSV for further analysis
	//
	// Included QC sources:
	//   - FastQC: Read quality metrics
	//   - SALMON: Quantification statistics
	//   - STAR: Alignment statistics
	//   - RSeQC: RNA-seq specific QC metrics
	// =====================================================================================
	
	multiqc_ch = Channel.empty()
		// FastQC outputs (ZIP files contain machine-readable data for MultiQC)
		// HTML files are for manual inspection only
		.mix(FASTQC_RAW.out.fastqc_zip)
		
		// SALMON outputs (only the directory path, not the sample_id)
		// .map { it[1] } extracts just the directory from [sample_id, directory] tuple
		// SALMON directory contains logs.dir/ with mapping statistics
		.mix(SALMON_QUANT.out.salmon_quant.map { it[1] })
		
		// STAR gene count files
		// Contains unstranded and stranded read counts per gene
		.mix(STAR_ALIGN.out.gene_counts)
		
		// STAR splice junction files
		// Contains all detected splice junctions
		.mix(STAR_ALIGN.out.sj_tab)
		
		// STAR alignment log files
		// Contains mapping rates and alignment statistics
		.mix(STAR_ALIGN.out.star_log)
		
		// RSeQC outputs (text, log, and R script files)
		// Contains read distribution, gene body coverage, etc.
		.mix(RSEQC.out.rseqc_logs)
		
		// Collect all items into a single list (waits for all samples)
		// This ensures MultiQC has access to all QC data before running
		.collect()
	
	// Run MultiQC with all collected QC files
	// Output will be a single HTML report in 06.MultiQC/ directory
	MULTIQC(multiqc_ch)
	
}

// =========================================================================================
// PIPELINE COMPLETION NOTES
// =========================================================================================
//
// 1. OUTPUT ORGANIZATION:
//    ${proj_dir}/
//    ├── 01.FastQ/
//    │   └── raw/                     # Raw FASTQ files (user-provided)
//    ├── 02.FastQC/
//    │   ├── raw/                     # FastQC reports on raw reads
//    │   │   ├── *.html               # Interactive HTML reports
//    │   │   └── *.zip                # Data files for MultiQC
//    │   └── trimmed/                 # FastQC reports on trimmed reads (if applicable)
//    ├── 03.Salmon/
//    │   ├── Sample1/                 # SALMON output directories
//    │   │   ├── quant.sf             # Main quantification file
//    │   │   ├── aux_info/            # Auxiliary information
//    │   │   └── logs/                # SALMON logs
//    │   ├── Sample2/
//    │   └── quant_files/             # Collected quant.sf files (all samples)
//    │       ├── Sample1.quant.sf
//    │       └── Sample2.quant.sf
//    ├── 04.STAR/
//    │   ├── Sample1.bam              # Aligned BAM files
//    │   ├── Sample1.bam.bai          # BAM indexes
//    │   ├── Sample1.ReadsPerGene.out.tab  # Gene counts
//    │   ├── Sample1.SJ.out.tab       # Splice junctions
//    │   └── Sample1.Log.final.out    # Alignment statistics
//    ├── 05.RSeQC/
//    │   ├── Sample1_read_distribution.txt        # Read distribution
//    │   ├── Sample1_geneBody_coverage.pdf        # Gene body coverage plot
//    │   ├── Sample1_junction_annotation.txt      # Junction annotation
//    │   └── Sample1_*.{pdf,txt,log,r}            # Other RSeQC outputs
//    ├── 06.MultiQC/
//    │   ├── Project_MultiQC_Report.html          # Final aggregated report
//    │   └── Project_MultiQC_Report_data/         # Data files (JSON, CSV)
//    │       ├── multiqc_data.json
//    │       └── multiqc_*.txt
//    └── 07.Logs/
//        ├── trace.txt                # Execution trace
//        ├── report.html              # Execution report
//        ├── timeline.html            # Timeline visualization
//        └── *.error.log              # Process error logs
//
// 2. DOWNSTREAM ANALYSIS:
//    a) Differential Expression with DESeq2/edgeR:
//       Option 1: Use STAR ReadsPerGene.out.tab files
//         - Load gene count matrices directly
//         - Choose appropriate column (unstranded, strand1, or strand2)
//         - Combine all samples into a count matrix
//       Option 2: Import SALMON quant.sf files with tximport
//         - Recommended for transcript-level analysis
//         - Handles multi-mapping reads better
//         - Provides bias-corrected estimates
//    
//    b) Visualization:
//       - Load BAM files into IGV for gene-level inspection
//         - Requires BAM and BAI files
//         - Can visualize splice junctions, coverage, etc.
//       - Generate coverage plots with RSeQC or deepTools
//         - Create heatmaps, profile plots
//         - Compare across samples
//    
//    c) Splice Variant Analysis:
//       - Use STAR SJ.out.tab files for junction analysis
//       - Consider tools like rMATS (differential splicing)
//       - Or LeafCutter (intron clustering)
//    
//    d) Gene Set Enrichment Analysis (GSEA):
//       - Use SALMON TPM values for ranking
//       - Or use DESeq2/edgeR results for ranked lists
//
// 3. QUALITY CONTROL CHECKPOINTS:
//    Review MultiQC report for:
//    - Mapping rate >70% (lower suggests contamination or wrong genome)
//    - Duplication rate <50% (higher suggests low library complexity)
//    - Consistent read distribution across samples
//      - Similar % in CDS, UTRs, introns across samples
//      - Batch effects show up as systematic differences
//    - Minimal 3' bias (ratio <3)
//      - High 3' bias suggests RNA degradation
//    - Low intergenic/intronic reads (<10%)
//      - High values suggest genomic DNA contamination
//    - Even gene body coverage
//      - Uniform coverage from 5' to 3' is ideal
//    - Consistent insert size (PE data)
//      - Should be consistent within expected range
//
// 4. RESUMING FAILED RUNS:
//    Nextflow caches completed processes in the work/ directory
//    To resume from last successful step:
//    
//    nextflow run main.nf -profile xinyi -resume
//    
//    This will:
//    - Skip successfully completed tasks
//    - Re-run only failed or incomplete tasks
//    - Use cached results from work/ directory
//    - Save significant computation time
//
// 5. CLEANING UP:
//    Remove work directory to free disk space after successful run:
//    
//    nextflow clean -f
//    
//    Or for more control:
//    nextflow clean -f -k              # Keep last run
//    nextflow clean -f -before <date>  # Clean runs before date
//    
//    Note: Only do this after confirming all outputs are correct!
//          You won't be able to resume the pipeline after cleaning.
//          Published outputs in project directories are not affected.
//
// 6. COMMON ISSUES AND SOLUTIONS:
//    a) "No such file or directory" errors:
//       - Check that raw_fastq_dir path is correct
//       - Verify reference genome files exist
//       - Ensure proper file permissions
//       - Check Singularity bind paths in nextflow.config
//    
//    b) Out of memory errors:
//       - Increase memory in process labels (nextflow.config)
//       - Check memory = { X.GB * task.attempt } settings
//       - Reduce parallelization with executor.queueSize
//       - Monitor with timeline.html to identify memory-hungry processes
//    
//    c) Process hangs indefinitely:
//       - Check system resource availability (disk space, memory)
//       - Review error logs in work/ directories
//       - Check .command.log files for specific errors
//       - Verify container images are accessible
//    
//    d) Container-related errors:
//       - Check Singularity cache directory exists
//       - Verify network access for container downloads
//       - Ensure bind paths match actual file locations
//       - Review runOptions in singularity configuration
//
// 7. EXTENDING THE PIPELINE:
//    To add new processes:
//    - Create new module file in modules/ directory
//      - Follow existing module structure
//      - Define inputs, outputs, script block
//      - Add appropriate labels for resource allocation
//    - Add process import at top of main.nf
//      - include { NEW_PROCESS } from './modules/new_process.nf'
//    - Insert process call in workflow block
//      - Place in appropriate position based on dependencies
//    - Add publishDir configuration in nextflow.config
//      - Define where outputs should be copied
//      - Set appropriate patterns and modes
//    - Update MultiQC input channel if process generates QC data
//      - Mix new QC outputs into multiqc_ch
//
// 8. PERFORMANCE OPTIMIZATION:
//    - Use -profile sge for cluster execution
//      - Parallelizes across multiple nodes
//      - Set appropriate queueSize to avoid overloading
//    - Monitor resource usage with timeline.html
//      - Identify bottlenecks
//      - Adjust CPU/memory allocations
//    - Consider trimming adapters if needed
//      - Add trimming process before alignment
//      - Can improve mapping rates
//    - Use -process.cache 'deep' for aggressive caching
//      - Useful during development
//      - Caches even if scripts change slightly
//
// =========================================================================================