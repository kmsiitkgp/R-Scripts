// =========================================================================================
// PROCESS: SALMON_QUANT
// =========================================================================================
// Purpose: Performs fast, alignment-free transcript quantification using SALMON
//
// SALMON advantages over traditional alignment:
//   - Much faster (~10-30x faster than STAR + featureCounts)
//   - Does not require alignment to genome
//   - Better handles multi-mapping reads through EM algorithm
//   - Provides transcript-level AND gene-level abundance estimates
//   - Corrects for biases: GC content, sequence-specific, positional
//
// Output files (per sample):
//   - quant.sf: Transcript abundances (TPM, NumReads)
//   - meta_info.json: Run metadata and parameters
//   - lib_format_counts.json: Library type detection results
//   - aux_info/: Auxiliary information for downstream tools
//
// Use cases:
//   - Differential gene expression (import to DESeq2/edgeR via tximport)
//   - Isoform-level analysis
//   - Quick QC before full alignment
// =========================================================================================

process SALMON_QUANT {
	
	tag "$sample_id"  // Display sample being quantified in logs
    label 'process_medium'  // Moderate CPU/memory requirements (4 cores, 12GB)
	
	// =================================================================================
	// INPUT
	// =================================================================================
	// Receives grouped sample tuples from VALIDATE_INPUT
	input:
	tuple val(sample_id), path(fastq_files)  // [sample_id, [R1.fq.gz] or [R1.fq.gz, R2.fq.gz]]
	path(salmon_index_dir)                   // Pre-built SALMON index from PREP_REFERENCE
		
	// =================================================================================
	// OUTPUT
	// =================================================================================
	// SALMON creates a directory named after the sample containing multiple files
	output:
    tuple val(sample_id), path("${sample_id}"),		emit: salmon_quant              // Full output directory
	path("*.{fa,txt}"),								emit: salmon_intermediate_files // Debugging files (optional)
    path("*.SALMON_QUANT.error.log"),				emit: salmon_error_log          // Process log
	
	// =================================================================================
	// SCRIPT SETUP
	// =================================================================================
	script:
	
	def LOG = "${sample_id}.SALMON_QUANT.error.log"
	
	// Determine sequencing mode and construct appropriate SALMON arguments
	// SALMON uses different flags for SE vs PE data:
	//   PE: --mates1 R1.fq.gz --mates2 R2.fq.gz
	//   SE: --unmatedReads R1.fq.gz
	def MATES_ARGS = fastq_files.size() == 2 ?
		"--mates1 ${fastq_files[0]} --mates2 ${fastq_files[1]}" :
		"--unmatedReads ${fastq_files[0]}"	
	
	"""	
	
	# =============================================================================
	# Run SALMON Quantification
	# =============================================================================
	# SALMON quantifies transcript abundances without alignment using k-mer matching
	#
	# Core parameters:
	#   --validateMappings: Use selective alignment mode (more accurate, slightly slower)
	#   --index: Path to SALMON index directory
	#   --output: Directory to store results (named after sample)
	#   --threads: Number of CPU cores to use
	#
	# Bias correction parameters (from params.SALMON_ARGS):
	#   --libType A: Auto-detect library type (strand orientation)
	#   --gcBias: Correct for GC content bias in fragments
	#   --seqBias: Correct for sequence-specific biases (e.g., random hexamer priming)
	#   --posBias: Correct for positional bias along transcripts
	#
	# These corrections improve quantification accuracy by 5-10%, especially for:
	#   - Genes with extreme GC content
	#   - Protocols using random hexamer priming
	#   - 3'-biased protocols
	# =============================================================================
	
	salmon quant \
		--validateMappings \
		--index "${salmon_index_dir}" \
		${MATES_ARGS} \
		${params.SALMON_ARGS.join(' ')} \
		--threads "${task.cpus}" \
		--output "${sample_id}" \
		1>> "${LOG}" 2>&1 \
		|| { echo "❌ ERROR: Salmon quantification failed for ${sample_id}" | tee -a "${LOG}" >&2; exit 1; }
		
	echo "✅ SUCCESS: Salmon quantification completed for ${sample_id}" >> "${LOG}"
		
	"""	
}

// =========================================================================================
// TECHNICAL NOTES
// =========================================================================================
//
// 1. SALMON QUANTIFICATION APPROACH:
//    - Uses k-mer based "quasi-mapping" instead of full alignment
//    - Builds a compatibility graph of reads to transcripts
//    - Uses Expectation-Maximization (EM) to assign multi-mapping reads
//    - Result: Probabilistic assignment of reads to transcripts
//
// 2. OUTPUT FILE: quant.sf
//    - TSV file with columns:
//      * Name: Transcript ID
//      * Length: Effective transcript length
//      * EffectiveLength: Length accounting for biases
//      * TPM: Transcripts Per Million (normalized abundance)
//      * NumReads: Estimated number of reads from this transcript
//    - This file is parsed by MultiQC for reporting
//
// 3. LIBRARY TYPE AUTO-DETECTION (--libType A):
//    - SALMON automatically determines strand orientation:
//      * SF: Single-end forward stranded
//      * SR: Single-end reverse stranded
//      * U: Unstranded
//      * ISF: Paired-end forward stranded (R2 matches transcript)
//      * ISR: Paired-end reverse stranded (R1 matches transcript)
//      * IU: Paired-end unstranded
//    - Detection results stored in lib_format_counts.json
//
// 4. BIAS CORRECTIONS EXPLAINED:
//    a) GC Bias (--gcBias):
//       - PCR amplification favors certain GC content
//       - SALMON models this bias and corrects abundances
//       - Most important for non-UMI protocols
//    
//    b) Sequence Bias (--seqBias):
//       - Random hexamer priming isn't truly random
//       - Creates bias toward specific 6-mer sequences
//       - SALMON learns and corrects for hexamer preferences
//    
//    c) Positional Bias (--posBias):
//       - Fragments aren't uniformly distributed along transcripts
//       - 3' bias in degraded samples or poly-A selection
//       - 5' bias in some library prep methods
//       - SALMON models coverage distribution and corrects
//
// 5. VALIDATEMAPPINGS MODE:
//    - More stringent than default "fast" mode
//    - Performs alignment of k-mer chains to validate mappings
//    - Reduces false positives from spurious k-mer matches
//    - Recommended for publication-quality results
//    - Trade-off: ~20% slower but more accurate
//
// 6. INTEGRATION WITH DOWNSTREAM ANALYSIS:
//    - DESeq2/edgeR: Use tximport package to import SALMON results
//      * tximport handles transcript->gene aggregation
//      * Preserves inferential uncertainty for better statistics
//    - Example R code:
//      ```r
//      library(tximport)
//      files <- list.files(pattern = "quant.sf", recursive = TRUE, full.names = TRUE)
//      names(files) <- basename(dirname(files))
//      txi <- tximport(files, type = "salmon", tx2gene = tx2gene)
//      dds <- DESeqDataSetFromTximport(txi, colData, ~condition)
//      ```
//
// 7. PUBLISHDIR CONFIGURATION (from nextflow.config):
//    - Full directory: Published to 03.Salmon/${sample_id}/
//    - quant.sf only: Copied to 03.Salmon/quant_files/ with renamed filename
//    - This allows easy collection of all quant.sf files for tximport
//
// 8. PERFORMANCE NOTES:
//    - Speed: ~2-5 minutes per sample on 4 cores
//    - Memory: ~8-12GB for human transcriptome
//    - Scales well with more cores (near-linear up to 8 cores)
//    - Bottleneck is usually I/O, not CPU
//
// 9. COMMON ISSUES AND SOLUTIONS:
//    - Low mapping rate (<50%):
//      * Wrong index (check species)
//      * Heavy contamination
//      * Degraded RNA (check FastQC)
//    - Very high duplication rate:
//      * Low RNA input
//      * Over-amplification
//      * Check library complexity
//    - Inconsistent library type:
//      * Mixed libraries in same run
//      * Check library prep protocol
// ========================================================================================='