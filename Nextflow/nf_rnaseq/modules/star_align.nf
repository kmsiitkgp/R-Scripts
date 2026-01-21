// =========================================================================================
// PROCESS: STAR_ALIGN
// =========================================================================================
// Purpose: Performs splice-aware alignment of RNA-seq reads to the reference genome
//          using the STAR (Spliced Transcripts Alignment to a Reference) aligner
//
// STAR Features:
//   - Ultra-fast: ~10-30 minutes per sample on 8 cores
//   - Splice-aware: Accurately maps reads across exon-exon junctions
//   - Two-pass mode: Discovers novel splice junctions, then re-aligns
//   - Outputs: Aligned BAM file + gene counts + splice junction table
//
// Key outputs:
//   - Aligned.sortedByCoord.out.bam: Coordinate-sorted BAM file for visualization/QC
//   - ReadsPerGene.out.tab: Gene-level read counts (for DESeq2/edgeR)
//   - SJ.out.tab: Splice junction coordinates and supporting reads
//   - Log.final.out: Alignment statistics (mapping rate, unique vs multi-mappers)
//
// Use cases:
//   - Gene-level differential expression (using ReadsPerGene.out.tab)
//   - Splice variant analysis (using SJ.out.tab)
//   - Visualization in IGV (using BAM file)
//   - QC analysis with RSeQC (using BAM file)
// =========================================================================================

process STAR_ALIGN {	
	
	tag "$sample_id"  // Display sample being aligned in logs
    label 'process_high'  // STAR alignment requires substantial RAM (30-50GB)
	
	// =================================================================================
	// INPUT
	// =================================================================================
	// Receives grouped sample tuples from VALIDATE_INPUT
	input:
	tuple val(sample_id), path(fastq_files)  // [sample_id, [R1.fq.gz] or [R1.fq.gz, R2.fq.gz]]
	path(star_index_dir)                     // Pre-built STAR index from PREP_REFERENCE
	
	// =================================================================================
	// OUTPUT
	// =================================================================================
	// STAR generates multiple output files with standardized naming
	output:
    tuple val(sample_id),
		path("${sample_id}.Aligned.sortedByCoord.out.bam"),      // Sorted BAM file
		path("${sample_id}.Aligned.sortedByCoord.out.bam.bai"),  // BAM index
		emit: bam_indexed  // Emit as tuple for downstream processes (e.g., RSeQC)

    path "${sample_id}.ReadsPerGene.out.tab",   emit: gene_counts     // Gene count matrix
    path "${sample_id}.SJ.out.tab",             emit: sj_tab          // Splice junction table
	path "${sample_id}.Log.final.out",          emit: star_log        // Alignment statistics
    path("*.STAR_ALIGN.error.log"),				emit: star_error_log  // Process log
		
	// =================================================================================
	// SCRIPT SETUP
	// =================================================================================
	script:	

	// Determine sequencing mode and construct appropriate STAR arguments
	// STAR uses different flags for SE vs PE data:
	//   PE: --readFilesIn R1.fq.gz R2.fq.gz
	//   SE: --readFilesIn R1.fq.gz
	def MATES_ARGS = fastq_files.size() == 2 ?
		"--readFilesIn ${fastq_files[0]} ${fastq_files[1]}" :
		"--readFilesIn ${fastq_files[0]}"	
	
	def LOG = "${sample_id}.STAR_ALIGN.error.log"
	
	"""	
	# =============================================================================
	# STEP 1: Align Reads Using STAR
	# =============================================================================
	# STAR performs splice-aware alignment in two passes:
	#   Pass 1: Align reads, discover novel splice junctions
	#   Pass 2: Re-align using both known (GTF) + novel junctions
	#
	# Core parameters:
	#   --genomeDir: Path to STAR index
	#   --readFilesIn: Input FASTQ file(s)
	#   --outFileNamePrefix: Prefix for all output files (sample_id.)
	#   --runThreadN: Number of CPU cores
	#
	# Additional parameters from params.STAR_ARGS (defined in nextflow.config):
	#   --runMode alignReads: Perform alignment (vs. index generation)
	#   --twopassMode Basic: Enable two-pass mode for novel junction discovery
	#   --quantMode GeneCounts: Generate gene count table
	#   --sjdbOverhang 100: Optimal for 101bp reads (ReadLength - 1)
	#   --readFilesCommand zcat: Decompress .gz files on-the-fly
	#   --outSAMtype BAM SortedByCoordinate: Output sorted BAM instead of SAM
	#   
	# See detailed parameter explanations in the notes section below
	# =============================================================================
	
	STAR \
		--genomeDir "${star_index_dir}" \
		${MATES_ARGS} \
		${params.STAR_ARGS.join(' ')} \
		--outFileNamePrefix "${sample_id}." \
		--runThreadN "${task.cpus}" \
		1>> "${LOG}" 2>&1 \
		|| { echo "❌ ERROR: STAR alignment failed for ${sample_id}" | tee -a "${LOG}" >&2; exit 1; }
	
	echo "✅ SUCCESS: STAR alignment completed for ${sample_id}" >> "${LOG}"
	
	# =============================================================================
	# STEP 2: Index the BAM File
	# =============================================================================
	# Generate BAI index file for fast random access to BAM
	# Required for:
	#   - IGV visualization
	#   - RSeQC analysis
	#   - Region-specific queries
	#
	# Using sambamba instead of samtools for speed (multi-threaded indexing)
	# =============================================================================
	
	sambamba index "${sample_id}.Aligned.sortedByCoord.out.bam" \
		1>> "${LOG}" 2>&1 \
		|| { echo "❌ ERROR: BAI index generation failed for ${sample_id}" | tee -a "${LOG}" >&2; exit 1; }
	
	echo "✅ SUCCESS: BAI index generation completed for ${sample_id}" >> "${LOG}"
	
	"""    
}

// =========================================================================================
// STAR COMMAND-LINE ARGUMENTS: DETAILED TECHNICAL EXPLANATIONS
// =========================================================================================
//
// This section explains each STAR parameter defined in nextflow.config (params.STAR_ARGS)
// Understanding these parameters is critical for optimizing alignment quality
//
// -----------------------------------------------------------------------------------------
// BASIC OPERATION PARAMETERS
// -----------------------------------------------------------------------------------------
//
// --runMode alignReads
//   Purpose: Tells STAR to perform read alignment (vs. genome index generation)
//   When to change: Never during alignment. Only use "genomeGenerate" when indexing
//
// --twopassMode Basic
//   Purpose: Enables two-pass alignment strategy
//   How it works:
//     Pass 1: STAR aligns reads and discovers novel splice junctions not in GTF
//     Pass 2: STAR rebuilds splice junction database + re-aligns all reads
//   Benefits:
//     - Discovers novel isoforms and fusion events
//     - Improves mapping rate by 2-5%
//     - Critical for organisms with incomplete annotations
//   Trade-off: ~30% slower than single-pass mode
//   Recommendation: Always use for RNA-seq unless speed is critical
//
// --quantMode GeneCounts
//   Purpose: Generates ReadsPerGene.out.tab with gene-level counts
//   Output columns:
//     Column 1: Gene ID
//     Column 2: Counts (unstranded)
//     Column 3: Counts (sense strand)
//     Column 4: Counts (antisense strand)
//   Use case: Direct input to DESeq2/edgeR for differential expression
//   Note: Also supports "TranscriptomeSAM" for isoform-level analysis
//
// --sjdbOverhang 100
//   Purpose: Sets the maximum overhang for splice junctions
//   Calculation: Ideally ReadLength - 1
//   Examples:
//     - 100bp reads → use 100
//     - 150bp reads → use 149
//     - 75bp reads → use 74
//   Impact: Affects splice junction detection sensitivity
//   Recommendation: Use 100 as default (works well for 75-150bp reads)
//
// -----------------------------------------------------------------------------------------
// FILE HANDLING PARAMETERS
// -----------------------------------------------------------------------------------------
//
// --readFilesCommand zcat
//   Purpose: Decompress .gz files during alignment
//   Alternatives:
//     - "bzcat" for .bz2 files
//     - "uncompress -c" for .Z files
//   Benefits: No need to pre-decompress FASTQ files (saves disk space)
//   Note: Slightly slower than pre-decompressed files (~5% overhead)
//
// --outSAMtype BAM SortedByCoordinate
//   Purpose: Output format and sorting
//   Options:
//     - SAM: Text format (large, slow)
//     - BAM: Binary compressed (small, fast)
//     - Unsorted: Fast output, but not usable for most tools
//     - SortedByCoordinate: Required for IGV, RSeQC, most downstream tools
//   Recommendation: Always use "BAM SortedByCoordinate"
//
// --outSAMunmapped Within
//   Purpose: Include unmapped reads in main BAM file
//   Alternatives:
//     - "None": Discard unmapped reads
//     - "Within KeepPairs": Keep pairs even if one mate unmapped
//   Benefits: Easier troubleshooting (can see why reads failed to map)
//   Trade-off: Slightly larger BAM files
//
// -----------------------------------------------------------------------------------------
// MULTI-MAPPING PARAMETERS
// -----------------------------------------------------------------------------------------
//
// --outFilterMultimapNmax 10
//   Purpose: Maximum number of loci a read can map to
//   Behavior: Reads mapping to >10 loci are discarded
//   Why allow multi-mappers?
//     - Gene families (e.g., HOX genes, immunoglobulins)
//     - Recent gene duplications
//     - Repetitive elements
//   Impact: Setting to 1 (unique-only) loses 10-20% of data
//   Recommendation: 10-20 for gene expression, 1 for variant calling
//
// -----------------------------------------------------------------------------------------
// MISMATCH FILTERING PARAMETERS
// -----------------------------------------------------------------------------------------
//
// --outFilterMismatchNmax 999
//   Purpose: Maximum total mismatches allowed per read
//   Why set so high: We use the percentage filter (below) instead
//   Strategy: Don't reject reads by total count, use percentage of read length
//
// --outFilterMismatchNoverReadLmax 0.04
//   Purpose: Maximum mismatches as fraction of read length
//   Calculation: 4% of read length
//   Examples:
//     - 100bp read: max 4 mismatches
//     - 150bp read: max 6 mismatches
//     - 250bp read: max 10 mismatches
//   Rationale: Longer reads naturally have more errors
//   Recommendation: 0.04-0.06 for most applications
//
// --outFilterType BySJout
//   Purpose: Filter reads based on splice junction confidence
//   How it works:
//     - Requires junctions to be either:
//       a) Present in GTF annotation, OR
//       b) Supported by multiple reads
//   Benefits: Reduces false novel junctions from sequencing errors
//   Trade-off: May miss genuine rare splice variants
//   Recommendation: Use for gene expression, consider disabling for isoform discovery
//
// -----------------------------------------------------------------------------------------
// ALIGNMENT STRATEGY PARAMETERS
// -----------------------------------------------------------------------------------------
//
// --alignEndsType Local
//   Purpose: Allow soft-clipping at read ends
//   How it works:
//     - STAR can "clip" (ignore) messy ends that don't align well
//     - Common causes: adapter contamination, low-quality bases
//   Benefits:
//     - Rescues reads that would otherwise fail mismatch filter
//     - Handles adapter read-through in short inserts
//   Alternative: "EndToEnd" (no clipping, more stringent)
//   Recommendation: Use "Local" for RNA-seq
//
// --alignIntronMin 20
//   Purpose: Minimum gap size to be considered an intron
//   Behavior: Gaps <20bp are classified as deletions (not introns)
//   Rationale:
//     - Biological introns are rarely <20bp
//     - Small gaps are likely sequencing errors or small deletions
//   Impact: Prevents spurious junction calls from deletions
//   Recommendation: 20-50bp for most organisms
//
// --alignIntronMax 1000000
//   Purpose: Maximum allowed intron size
//   Behavior: Gaps >1Mb are not considered valid introns
//   Rationale:
//     - Human/mouse largest introns: ~500-800kb
//     - Prevents misalignment across chromosomes
//   Recommendation: 1000000 (1Mb) for mammals
//
// --alignMatesGapMax 1000000
//   Purpose: Maximum distance between paired-end mates
//   Behavior: Read pairs >1Mb apart are considered discordant
//   Use case: Should match --alignIntronMax for RNA-seq
//   Note: For DNA-seq, would be much smaller (typically 1-2kb)
//
// -----------------------------------------------------------------------------------------
// PERFORMANCE CONSIDERATIONS
// -----------------------------------------------------------------------------------------
//
// Resource requirements:
//   - Memory: 30-50GB for human genome (index loaded into RAM)
//   - CPU: Scales well up to 8-16 cores
//   - Time: ~10-30 min/sample on 8 cores
//   - Disk: ~5-10GB per sample for output BAM
//
// Memory optimization:
//   - Use --limitBAMsortRAM to control BAM sorting memory
//   - Use --outBAMsortingThreadN for parallel sorting
//   - For low memory: use --genomeLoad NoSharedMemory
//
// Speed optimization:
//   - More cores = faster (diminishing returns >16 cores)
//   - SSD storage significantly faster than HDD
//   - Pre-decompressed FASTQ slightly faster than .gz
//
// -----------------------------------------------------------------------------------------
// OUTPUT FILES EXPLAINED
// -----------------------------------------------------------------------------------------
//
// Aligned.sortedByCoord.out.bam
//   - Main output: coordinate-sorted BAM file
//   - Contains: aligned reads with CIGAR strings, flags, quality scores
//   - Use for: visualization (IGV), QC (RSeQC), variant calling
//   - Size: 2-10GB per sample
//
// ReadsPerGene.out.tab
//   - Gene-level read counts (4 columns)
//   - Column 2 (unstranded): Use for unstranded libraries
//   - Column 3 (sense): Use for stranded libraries (dUTP, Illumina)
//   - Column 4 (antisense): Use for reverse-stranded protocols
//   - Direct input to DESeq2/edgeR
//
// SJ.out.tab
//   - Splice junction table (9 columns):
//     1. Chromosome
//     2. Start position (intron)
//     3. End position (intron)
//     4. Strand (0=undefined, 1=+, 2=-)
//     5. Motif (0=non-canonical, 1=GT/AG, 2=CT/AC, etc.)
//     6. Annotation (0=unannotated, 1=annotated)
//     7. Uniquely mapping reads crossing junction
//     8. Multi-mapping reads crossing junction
//     9. Maximum overhang
//   - Use for: splice variant analysis, fusion detection
//
// Log.final.out
//   - Alignment statistics:
//     * Total reads
//     * Uniquely mapped reads (%)
//     * Reads mapped to multiple loci (%)
//     * Unmapped reads (%)
//     * Mismatch rate
//     * Deletion/insertion rates
//   - Critical for QC: aim for >70% unique mapping
//
// -----------------------------------------------------------------------------------------
// TROUBLESHOOTING COMMON ISSUES
// -----------------------------------------------------------------------------------------
//
// Low mapping rate (<50%):
//   - Wrong genome/index (check species)
//   - Heavy contamination (bacteria, adapter)
//   - Degraded RNA (check RIN scores)
//   - Wrong library type (DNA instead of RNA)
//
// High multi-mapping rate (>30%):
//   - Ribosomal RNA contamination
//   - Low library complexity
//   - Repetitive sequences
//   - Consider rRNA depletion
//
// Out of memory errors:
//   - Increase memory allocation
//   - Use --limitBAMsortRAM to cap memory
//   - Use smaller --genomeSAindexNbases during indexing
//
// Slow performance:
//   - Increase CPU cores
//   - Use SSD storage
//   - Check I/O bottlenecks
//   - Consider disabling two-pass mode
// ========================================================================================='