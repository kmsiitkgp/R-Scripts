// =========================================================================================
// PROCESS: RSEQC
// =========================================================================================
// Purpose: Performs comprehensive quality control analysis on aligned RNA-seq BAM files
//          using the RSeQC package
//
// RSeQC is a collection of Python scripts that provide detailed QC metrics including:
//   - Read distribution across genomic features (CDS, UTR, introns, intergenic)
//   - Insert size distribution (paired-end only)
//   - Splice junction annotation (known vs novel)
//   - Junction saturation analysis
//   - Gene body coverage uniformity (5' to 3' bias detection)
//   - Sequencing artifact profiling (deletions, mismatches, insertions, clipping)
//
// These metrics help identify:
//   - Library prep issues (3' bias, degradation)
//   - Contamination (non-coding RNA, genomic DNA)
//   - Sequencing quality problems
//   - Strand-specificity verification
// =========================================================================================

process RSEQC {
	
	tag "$sample_id"  // Display sample being analyzed in logs
    label 'process_medium'  // RSeQC can be memory-intensive on large BAMs (4 cores, 12GB)
	
	// =================================================================================
	// INPUT
	// =================================================================================
	input:	
    tuple val(sample_id), path(bam), path(bai)  // [sample_id, BAM_file, BAI_index]
	path(ref_bed)                                // Gene annotation in BED format
	val(mode)                                    // "SINGLE_END" or "PAIRED_END"
	
	// =================================================================================
	// OUTPUT
	// =================================================================================
	// RSeQC generates multiple plot and data files
	// Some outputs are only generated for paired-end data (e.g., inner_distance)
	output:
	path("${sample_id}*.{pdf,jpeg,png,tiff}"),  emit: rseqc_plots,       optional: true  // Plots
    path("${sample_id}*.{txt,log,r}"),          emit: rseqc_logs,        optional: true  // Data files
    path("*.RSEQC.error.log"),					emit: rseqc_error_log                    // Process log

    // =================================================================================
	// SCRIPT SETUP
	// =================================================================================
	script:
	
    // Translate Nextflow 'mode' to RSeQC 'sequencing' parameter
    // RSeQC uses "PE" or "SE" for some scripts
    def SEQUENCING_MODE = (mode == "PAIRED_END") ? "PE" : "SE"

	def LOG = "${sample_id}.RSEQC.error.log"
	
	"""
	# =============================================================================
	# Calculate Read Length from BAM
	# =============================================================================
	# Some RSeQC scripts require the read length as input
	# Extract from first read in BAM using samtools + awk
	# Column 10 in SAM format contains the sequence
	# =============================================================================
	
	READ_LENGTH=\$(samtools view "${bam}" | head -n1 | awk '{print length(\$10)}')
	
	# =============================================================================
	# QC STEP 1: Calculate Read Distribution Across Genomic Features
	# =============================================================================
	# Script: read_distribution.py
	# Purpose: Determines where reads map in the genome
	# 
	# Output categories:
	#   - CDS (coding sequence): Expected for most mRNA
	#   - 5' UTR: Should be present
	#   - 3' UTR: Should be present (may be enriched in 3' biased protocols)
	#   - Introns: Should be low (<10%) for good RNA quality
	#   - Intergenic: Should be low (<5%)
	#
	# Interpretation:
	#   - High intron %: Genomic DNA contamination or poor RNA quality
	#   - High intergenic %: Wrong genome annotation or contamination
	#   - Low CDS %: Degraded RNA or strand-specificity issues
	# =============================================================================
	
	read_distribution.py \
		--input-file "${bam}" \
		--refgene "${ref_bed}" \
		> "${sample_id}.RSeQC.txt" \
		1>> "${LOG}" 2>&1 \
		|| { echo "❌ ERROR: Read distribution calculation failed for ${sample_id}" | tee -a "${LOG}" >&2; exit 1; }
	
	echo "✅ SUCCESS: Read distribution calculation completed for ${sample_id}" >> "${LOG}"

	
	# =============================================================================
	# QC STEP 2: Calculate Inner Distance Between Read Pairs (PAIRED-END ONLY)
	# =============================================================================
	# Script: inner_distance.py
	# Purpose: Measures the distance between paired-end mates (insert size - read length)
	# 
	# Expected range: 50-300bp for typical mRNA libraries
	# 
	# Outputs:
	#   - *.inner_distance.txt: Distance statistics
	#   - *.inner_distance_plot.pdf: Distribution histogram
	#
	# Interpretation:
	#   - Very large distances: May indicate genomic DNA contamination
	#   - Very small distances: Adapter dimers or short inserts
	#   - Bimodal distribution: Mixed library or size selection issues
	# =============================================================================
	
	if [[ "${mode}" == "PAIRED_END" ]]; then
		inner_distance.py \
			--input-file "${bam}" \
			--refgene "${ref_bed}" \
			--mapq 30 \
			--out-prefix "${sample_id}" \
			1>> "${LOG}" 2>&1 \
			|| { echo "❌ ERROR: Inner distance calculation failed for ${sample_id}" | tee -a "${LOG}" >&2; exit 1; }
		
		echo "✅ SUCCESS: Inner distance calculation completed for ${sample_id}" >> "${LOG}"
	fi
	
	# =============================================================================
	# QC STEP 3: Annotate Splice Junctions (Compare to Reference)
	# =============================================================================
	# Script: junction_annotation.py
	# Purpose: Classifies splice junctions as known (annotated) vs novel
	# 
	# Outputs:
	#   - *.junction.txt: Junction statistics
	#   - *.splice_junction.pdf: Pie chart of junction categories
	#   - *.splice_events.pdf: Event types (exon skipping, etc.)
	#
	# Categories:
	#   - Annotated: Both splice sites in GTF
	#   - Partial novel: One splice site in GTF
	#   - Complete novel: Neither splice site in GTF
	#
	# Interpretation:
	#   - >80% annotated: Good quality, expected for well-studied organisms
	#   - High novel %: Novel isoforms OR wrong annotation OR contamination
	#   - Very low junction count: Genomic DNA contamination
	# =============================================================================
	
	junction_annotation.py \
		--input-file "${bam}" \
		--refgene "${ref_bed}" \
		--min-intron 50 \
		--mapq 30 \
		--out-prefix "${sample_id}" \
		1>> "${LOG}" 2>&1 \
		|| { echo "❌ ERROR: Junction annotation failed for ${sample_id}" | tee -a "${LOG}" >&2; exit 1; }
	
	echo "✅ SUCCESS: Junction annotation completed for ${sample_id}" >> "${LOG}"		
	
	# =============================================================================
	# QC STEP 4: Check Junction Detection Saturation
	# =============================================================================
	# Script: junction_saturation.py
	# Purpose: Determines if sequencing depth is sufficient to detect junctions
	# 
	# How it works:
	#   - Subsamples reads at increasing depths (5%, 10%, ..., 100%)
	#   - Counts junctions detected at each depth
	#   - Plots saturation curve
	#
	# Outputs:
	#   - *.junctionSaturation_plot.pdf: Saturation curve
	#
	# Interpretation:
	#   - Plateau reached: Sufficient sequencing depth
	#   - Linear curve: Under-sequenced, more depth would detect more junctions
	#   - Recommended: >20M reads per sample for human
	# =============================================================================
	
	junction_saturation.py \
		--input-file "${bam}" \
		--refgene "${ref_bed}" \
		--min-intron 50 \
		--mapq 30 \
		--out-prefix "${sample_id}" \
		1>> "${LOG}" 2>&1 \
		|| { echo "❌ ERROR: Junction saturation calculation failed for ${sample_id}" | tee -a "${LOG}" >&2; exit 1; }
	
	echo "✅ SUCCESS: Junction saturation calculation completed for ${sample_id}" >> "${LOG}"

		
	# =============================================================================
	# QC STEP 5: Profile Deletion Rates by Read Position
	# =============================================================================
	# Script: deletion_profile.py
	# Purpose: Checks if deletions are uniformly distributed along reads
	# 
	# Outputs:
	#   - *.deletion_profile.txt: Deletion counts per position
	#   - *.deletion_profile.pdf: Line plot of deletion rate
	#
	# Interpretation:
	#   - Uniform distribution: Normal
	#   - Enrichment at ends: Sequencing quality degradation
	#   - Spikes at specific positions: Systematic sequencing artifacts
	# =============================================================================
	
	deletion_profile.py \
		--input "${bam}" \
		--read-align-length \${READ_LENGTH} \
		--mapq 30 \
		--out-prefix "${sample_id}" \
		1>> "${LOG}" 2>&1 \
		|| { echo "❌ ERROR: Deletion profile calculation failed for ${sample_id}" | tee -a "${LOG}" >&2; exit 1; }
	
	echo "✅ SUCCESS: Deletion profile calculation completed for ${sample_id}" >> "${LOG}"
		
	
	# =============================================================================
	# QC STEP 6: Profile Mismatch Rates by Read Position
	# =============================================================================
	# Script: mismatch_profile.py
	# Purpose: Checks if mismatches (SNPs/errors) are uniformly distributed
	# 
	# Outputs:
	#   - *.mismatch_profile.txt: Mismatch counts per position
	#   - *.mismatch_profile.pdf: Line plot of mismatch rate
	#
	# Interpretation:
	#   - Increasing toward 3' end: Quality degradation (normal for Illumina)
	#   - High mismatch rate overall: Wrong reference or contamination
	#   - Spikes at specific positions: Systematic errors or RNA editing
	# =============================================================================
	
	mismatch_profile.py \
		--input "${bam}" \
		--read-align-length \${READ_LENGTH} \
		--mapq 30 \
		--out-prefix "${sample_id}" \
		1>> "${LOG}" 2>&1 \
		|| { echo "❌ ERROR: Mismatch profile calculation failed for ${sample_id}" | tee -a "${LOG}" >&2; exit 1; }
	
	echo "✅ SUCCESS: Mismatch profile calculation completed for ${sample_id}" >> "${LOG}"
	
	
	# =============================================================================
	# QC STEP 7: Profile Insertion Rates by Read Position
	# =============================================================================
	# Script: insertion_profile.py
	# Purpose: Checks if insertions are uniformly distributed along reads
	# 
	# Outputs:
	#   - *.insertion_profile.txt: Insertion counts per position
	#   - *.insertion_profile.pdf: Line plot of insertion rate
	#
	# Interpretation:
	#   - Low, uniform distribution: Normal
	#   - High insertion rate: Potential contamination or wrong reference
	#   - Spikes: Systematic sequencing artifacts
	# =============================================================================
	
	insertion_profile.py \
		--input "${bam}" \
		--sequencing ${SEQUENCING_MODE} \
		--mapq 30 \
		--out-prefix "${sample_id}" \
		1>> "${LOG}" 2>&1 \
		|| { echo "❌ ERROR: Insertion profile calculation failed for ${sample_id}" | tee -a "${LOG}" >&2; exit 1; }
	
	echo "✅ SUCCESS: Insertion profile calculation completed for ${sample_id}" >> "${LOG}"
		
	
	# =============================================================================
	# QC STEP 8: Profile Soft-Clipping Rates at Read Ends
	# =============================================================================
	# Script: clipping_profile.py
	# Purpose: Checks how often read ends are soft-clipped (not aligned)
	# 
	# Outputs:
	#   - *.clipping_profile.txt: Clipping statistics
	#   - *.clipping_profile.pdf: Bar plots of clipping rates
	#
	# Interpretation:
	#   - Low clipping (<5%): Good quality
	#   - High clipping at 3' end: Adapter contamination
	#   - High clipping overall: Poor quality or complex library
	# =============================================================================
	
	clipping_profile.py \
		--input "${bam}" \
		--sequencing ${SEQUENCING_MODE} \
		--mapq 30 \
		--out-prefix "${sample_id}" \
		1>> "${LOG}" 2>&1 \
		|| { echo "❌ ERROR: Clipping profile calculation failed for ${sample_id}" | tee -a "${LOG}" >&2; exit 1; }
	
	echo "✅ SUCCESS: Clipping profile calculation completed for ${sample_id}" >> "${LOG}"
		
	# =============================================================================
	# QC STEP 9: Calculate Gene Body Coverage Uniformity (5' to 3' Bias)
	# =============================================================================
	# Script: geneBody_coverage.py
	# Purpose: Checks if coverage is uniform across gene bodies (5' to 3')
	# 
	# Process:
	#   1. Subsample BAM to 1M reads for speed (using samtools -s)
	#   2. Calculate coverage in 100 bins along gene bodies
	#   3. Average across all genes
	#
	# Outputs:
	#   - *.geneBodyCoverage.txt: Coverage values for each percentile
	#   - *.geneBodyCoverage.pdf: Line plot showing coverage profile
	#
	# Interpretation:
	#   - Flat line (uniform coverage): Ideal
	#   - 3' bias (increasing toward 3' end): Common in degraded RNA or poly-A protocols
	#   - 5' bias (increasing toward 5' end): Rare, may indicate specific library prep
	#   - Recommended: 3'/5' ratio <2 for good quality
	#
	# Note: 
	#   - Subsampling to 1M reads significantly speeds up analysis
	#   - -s 1.05: Random seed (1) + fraction (0.05 = 5%, but we want more so use 1.0 = 100%)
	#   - For 1M reads from 50M total: -s 1.02 (2% of reads)
	# =============================================================================
	
	# Subsample BAM to reduce processing time
	# -s seed.fraction: Random seed (1.0) + fraction needed for ~1M reads
	# For typical RNA-seq (20-50M reads), 5% gives us 1-2.5M reads
	samtools view -b -s 1.05 "${bam}" > "${sample_id}.subset.bam"
	
	geneBody_coverage.py \
		--input "${sample_id}.subset.bam" \
		--refgene "${ref_bed}" \
		--minimum_length 100 \
		--format pdf \
		--out-prefix "${sample_id}" \
		1>> "${LOG}" 2>&1 \
		|| { echo "❌ ERROR: Gene body coverage calculation failed for ${sample_id}" | tee -a "${LOG}" >&2; exit 1; }
	
	# Clean up temporary subset BAM
	rm "${sample_id}.subset.bam"
	
	echo "✅ SUCCESS: Gene body coverage calculation completed for ${sample_id}" >> "${LOG}"	

	"""
}

// =========================================================================================
// TECHNICAL NOTES
// =========================================================================================
//
// 1. WHY RSEQC IS CRITICAL:
//    - FastQC only checks raw read quality
//    - RSeQC checks if alignment makes biological sense
//    - Catches issues that appear only after alignment
//
// 2. COMMON QC ISSUES DETECTED:
//    a) High intron reads (>20%):
//       - Cause: Genomic DNA contamination or poor RNA extraction
//       - Fix: Improve RNA purification, use DNase treatment
//    
//    b) 3' bias (coverage ratio >3):
//       - Cause: RNA degradation or poly-A selection protocol
//       - Impact: May miss 5' transcript variants
//       - Fix: Improve RNA handling, use rRNA depletion instead
//    
//    c) Low junction detection:
//       - Cause: Genomic DNA contamination
//       - Fix: DNase treatment or better RNA extraction
//    
//    d) High novel junction rate (>30%):
//       - Cause: Wrong GTF annotation or contamination
//       - Fix: Verify correct genome/annotation version
//
// 3. PERFORMANCE OPTIMIZATION:
//    - Gene body coverage is the slowest step (~10-30 min)
//    - Subsampling to 1M reads reduces time by 10-50x
//    - Maintains statistical power for bias detection
//    - Can adjust subsampling rate based on sequencing depth
//
// 4. OUTPUT INTERPRETATION:
//    - RSeQC outputs are parsed by MultiQC for aggregation
//    - Individual PDFs useful for sample-specific investigation
//    - Text files contain raw data for custom analysis
//
// 5. PAIRED-END VS SINGLE-END:
//    - Only PE data gets inner_distance analysis
//    - Some scripts require --sequencing parameter (PE/SE)
//    - SE data runs faster (fewer reads to process)
//
// 6. MAPQ FILTERING (--mapq 30):
//    - Only uses high-quality alignments
//    - MAPQ 30 = 99.9% confidence in alignment
//    - Prevents contamination from multi-mappers
//
// 7. MINIMUM INTRON SIZE (--min-intron 50):
//    - Gaps <50bp considered deletions, not introns
//    - Matches biological reality (most introns >50bp)
//    - Prevents artifact junctions from indels
// ========================================================================================='