// =========================================================================================
// PROCESS: PREP_REFERENCE
// =========================================================================================
// Purpose: Prepares all necessary reference genome indexes for RNA-seq alignment and
//          quantification. This process is computationally expensive but only needs to
//          run once per genome version.
//
// Creates three key outputs:
//   1. STAR index - for splice-aware genome alignment
//   2. SALMON index - for fast transcript quantification
//   3. BED file - for RSeQC quality control analysis
//
// Note: These indexes are species-specific and published to the reference directory
//       for reuse across multiple projects
// =========================================================================================

process PREP_REFERENCE {
	
	tag "Indexing ${ref_fasta.baseName}"  // Display genome being indexed
    label 'process_high'  // STAR indexing requires substantial RAM (30-50GB)
	
	// =================================================================================
	// INPUT
	// =================================================================================
	// Reference genome files (species-specific, from nextflow.config)
	input:
	path(ref_fasta)  // Genome FASTA file (e.g., Homo_sapiens.GRCh38.dna.primary_assembly.fa)
    path(ref_gtf)    // Gene annotation GTF file (e.g., Homo_sapiens.GRCh38.113.gtf)
	
	// =================================================================================
	// OUTPUT
	// =================================================================================
	// All outputs are directories/files published to the reference genome directory
	// This allows reuse across multiple pipeline runs without re-indexing
	output:
	path("star_index_dir"),     		emit: star_index_dir             // STAR genome index
	path("salmon_index_dir"),   		emit: salmon_index_dir           // SALMON transcript index
	path("${ref_bed}"),         		emit: ref_bed                    // BED file for RSeQC
	path("PREP_REFERENCE.error.log"),	emit: prep_reference_error_log   // Process log
	
	// =================================================================================
	// SCRIPT SETUP
	// =================================================================================
	script:
	
	// Extract base name from GTF file for BED file naming
	// Example: "Homo_sapiens.GRCh38.113.gtf" -> "Homo_sapiens.GRCh38.113"
	def gtf_basename = ref_gtf.baseName 
	def ref_bed = "${gtf_basename}.bed"	
	
	def LOG = "PREP_REFERENCE.error.log"
	
	"""
	# =============================================================================
	# STEP 1: Create STAR Genome Index
	# =============================================================================
	# STAR requires a pre-built index for fast splice-aware alignment
	# 
	# Key parameters:
	#   --runMode genomeGenerate: Tells STAR to build an index (not align)
	#   --runThreadN: Number of parallel threads (from task.cpus)
	#   --genomeDir: Output directory for index files
	#   --genomeFastaFiles: Reference genome sequence
	#   --sjdbGTFfile: Gene annotations for splice junction database
	#   --sjdbOverhang: ReadLength-1 (100 is optimal for 101bp reads)
	#   --genomeSAindexNbases: Log2(genome_size)/2 - 1 (14 is optimal for human/mouse)
	#
	# Memory requirement: ~30GB for human genome
	# Time: ~30-60 minutes on 8 cores
	# =============================================================================
	
	mkdir -p star_index_dir
	
	STAR --runMode genomeGenerate \
        --runThreadN "${task.cpus}" \
        --genomeDir star_index_dir \
        --genomeFastaFiles "${ref_fasta}" \
        --sjdbGTFfile "${ref_gtf}" \
        --sjdbOverhang 100 \
        --genomeSAindexNbases 14 \
		1>> "${LOG}" 2>&1 \
		|| { echo "❌ ERROR: STAR index generation failed" | tee -a "${LOG}" >&2; exit 1; }
	
	echo "✅ SUCCESS: STAR index generation completed" >> "${LOG}"
		
	# =============================================================================
	# STEP 2: Create BED File for RSeQC
	# =============================================================================
	# RSeQC requires gene annotations in BED format (not GTF)
	# 
	# Conversion process:
	#   GTF -> genePred (intermediate format) -> BED (final format)
	#
	# Tools used:
	#   gtfToGenePred: UCSC tool to convert GTF to genePred format
	#   genePredToBed: UCSC tool to convert genePred to BED format
	#
	# The BED file contains:
	#   - Chromosome, start, end positions
	#   - Gene names and strand information
	#   - Exon/intron boundaries
	# =============================================================================
	
	gtfToGenePred "${ref_gtf}" tmp_file
	genePredToBed tmp_file "${ref_bed}" \
		1>> "${LOG}" 2>&1 \
		|| { echo "❌ ERROR: RSEQC BED generation failed" | tee -a "${LOG}" >&2; exit 1; }
	
	echo "✅ SUCCESS: RSEQC BED generation completed" >> "${LOG}"

	# =============================================================================
	# STEP 3: Create SALMON Index
	# =============================================================================
	# SALMON performs fast, alignment-free transcript quantification using k-mers
	# 
	# Process overview:
	#   1. Create decoy list from genome chromosomes
	#   2. Extract transcript sequences from genome using GTF annotations
	#   3. Combine transcripts + genome into "gentrome" (transcriptome + genome)
	#   4. Build k-mer index with decoy-aware mode
	#
	# Why use decoys?
	#   - Prevents false mapping of reads to transcripts when they actually 
	#     originate from intergenic/intronic regions
	#   - Improves quantification accuracy by ~10-15%
	#
	# Key parameters:
	#   --transcripts: Combined transcriptome and genome sequences (gentrome)
	#   --decoys: List of chromosome names to treat as decoy sequences
	#   --kmerLen 31: K-mer length for indexing (default, works well for reads >75bp)
	#   --threads: Number of parallel threads
	# =============================================================================
	
	# 3a. Extract chromosome names from FASTA headers to create decoy list
	# Example header: ">1 dna:chromosome chromosome:GRCh38:1:1:248956422:1"
	# Extracts: "1" (chromosome name)
	grep "^>" "${ref_fasta}" | cut -d " " -f1 | sed 's/>//' > salmon_decoy
	
	# 3b. Extract transcript sequences from genome using GTF coordinates
	# gffread reads GTF annotations and extracts corresponding sequences from FASTA
	gffread "${ref_gtf}" \
		-g "${ref_fasta}" \
		-w salmon_transcriptome_fasta	

	# 3c. Create gentrome by concatenating transcriptome + genome
	# IMPORTANT: Transcripts must come FIRST, genome sequences SECOND
	cat salmon_transcriptome_fasta "${ref_fasta}" > salmon_gentrome_fasta
	
	# 3d. Build the SALMON index with decoy-aware mode
	salmon index \
        --transcripts salmon_gentrome_fasta \
		--decoys salmon_decoy \
		--threads "${task.cpus}" \
		--kmerLen 31 \
		--index salmon_index_dir \
		1>> "${LOG}" 2>&1 \
		|| { echo "❌ ERROR: SALMON index generation failed" | tee -a "${LOG}" >&2; exit 1; }
	
	echo "✅ SUCCESS: SALMON index generation completed" >> "${LOG}"	
	
	"""	
}

// =========================================================================================
// TECHNICAL NOTES
// =========================================================================================
//
// 1. WHEN DOES THIS PROCESS RUN?
//    - First time you run the pipeline with a new genome version
//    - Index files are published to ref_dir and reused on subsequent runs
//    - Nextflow caching prevents re-running if inputs haven't changed
//
// 2. RESOURCE REQUIREMENTS:
//    - Memory: 30-50GB for human/mouse genomes
//    - CPU: 8+ cores recommended (indexing is parallelizable)
//    - Time: 1-2 hours for full human genome
//    - Disk: ~30GB for STAR index, ~5GB for SALMON index
//
// 3. STAR INDEX PARAMETERS EXPLAINED:
//    - sjdbOverhang: Should be ReadLength-1 (100 for 101bp reads)
//      * Affects splice junction detection sensitivity
//      * Using 100 works well for most read lengths (75-150bp)
//    - genomeSAindexNbases: Controls memory/speed trade-off
//      * Formula: min(14, log2(GenomeSize)/2 - 1)
//      * For human (3Gb): log2(3*10^9)/2 - 1 ≈ 14
//      * Smaller genomes may need smaller values (e.g., 13 for mouse)
//
// 4. SALMON GENTROME STRATEGY:
//    - Why combine transcriptome + genome?
//      * Allows SALMON to identify reads from intergenic regions
//      * Prevents false quantification of unannotated transcripts
//    - Decoy sequences (chromosomes) absorb reads that don't map to transcripts
//    - This "selective alignment" mode is now SALMON's recommended approach
//
// 5. BED FILE USAGE:
//    - Used by RSeQC scripts for QC analysis:
//      * read_distribution.py - where do reads map (CDS, UTR, intron, intergenic)?
//      * junction_annotation.py - are splice junctions known or novel?
//      * geneBody_coverage.py - is there 3' bias in coverage?
//    - BED format is simpler than GTF, making RSeQC scripts faster
//
// 6. ERROR HANDLING:
//    - Each major step has error checking with descriptive messages
//    - Logs capture both stdout and stderr for debugging
//    - Process fails immediately if any step fails (fail-fast approach)
//
// 7. PUBLISHDIR CONFIGURATION:
//    - Defined in nextflow.config
//    - Outputs go to: ${ref_dir}/${species}/
//    - Example: ~/NGSTools/Reference_Genomes/Human/star_index_dir
// ========================================================================================='