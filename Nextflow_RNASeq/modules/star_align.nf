process STAR_ALIGN {	
	
	def LOG = "${sample_id}.STAR_ALIGN.error.log"
	
	input:
	tuple val(sample_id), path(fastq_files)
	path(star_index_dir)
	
	output:
    tuple val(sample_id),
		path("${sample_id}.Aligned.sortedByCoord.out.bam"),
		path("${sample_id}.Aligned.sortedByCoord.out.bam.bai"),
		emit: bam_indexed 

    path "${sample_id}.ReadsPerGene.out.tab",		emit: gene_counts
    path "${sample_id}.SJ.out.tab",					emit: sj_tab    
	path "${sample_id}.Log.final.out",				emit: star_log
    path "${LOG}",							 		emit: star_error_log
		
	script:	

	// Determine if paired-end or single-end (Groovy, outside bash)
	def MATES_ARGS = fastq_files.size() == 2 ?
		"--readFilesIn ${fastq_files[0]} ${fastq_files[1]}" :
		"--readFilesIn ${fastq_files[0]}"	
	
	"""	
	
	# Align reads using STAR
	STAR \
		--genomeDir "${star_index_dir}" \
		${MATES_ARGS} \
		${params.STAR_ARGS.join(' ')} \
		--outFileNamePrefix "${sample_id}." \
		--runThreadN "${task.cpus}" \
		1>> "${LOG}" 2>&1 \
		|| { echo "❌ ERROR: STAR alignment failed for ${sample_id}" | tee -a "${LOG}" >&2; exit 1; }
	
	echo "✅ SUCCESS: STAR alignment completed for ${sample_id}" >> "${LOG}"
	
	# Generate index for BAM
	sambamba index "${sample_id}.Aligned.sortedByCoord.out.bam" \
		1>> "${LOG}" 2>&1 \
		|| { echo "❌ ERROR: BAI index generation failed for ${sample_id}" | tee -a "${LOG}" >&2; exit 1; }
	
	echo "✅ SUCCESS: BAI index generation completed for ${sample_id}" >> "${LOG}"
	
	"""    
}

/*	
# ==============================================================================
# 🚀 STAR COMMAND-LINE ARGUMENTS: TECHNICAL DEFINITIONS
# ==============================================================================

# --runMode alignReads
# Tells STAR to perform the actual alignment of input FASTQ reads to the 
# reference genome index.

# --twopassMode Basic
# Enables 2-pass mapping where STAR discovers novel splice junctions in the first 
# pass and uses them to guide a more accurate re-alignment in the second pass.
# Benefit: Critical for RSEM to correctly assign reads to novel isoforms.


# --quantMode TranscriptomeSAM GeneCounts
# 1. TranscriptomeSAM: Generates a BAM file in transcript coordinates. STAR 
#    automatically forces 'End-to-End' logic here to satisfy RSEM requirements.
# 2. GeneCounts: Generates 'ReadsPerGene.out.tab', counting overlaps with 
#    genomic features. Used for standard DGE (DESeq2/edgeR).

# --sjdbOverhang 100
# This value should ideally be (ReadLength - 1). For 101bp reads, 100 is optimal. 
# It provides the aligner with the correct "window" size to search for 
# junctions during the second pass of 2-pass mode.

# --readFilesCommand zcat
# Tells STAR to use zcat to open .gz files, allowing the use of compressed 
# FASTQ files directly without manual unzipping.

# --outFilterMultimapNmax 10
# Allows a read to map to up to 10 loci. If it maps to >10, it is discarded.
# Q: Why not unique-only (1)? 
# A: RSEM requires multi-mappers to statistically distribute counts among 
#    isoforms/paralogs, preventing significant data loss in gene-heavy regions.

# --outFilterMismatchNmax 999
# By setting this very high, we ensure that the total number of mismatches doesn't 
# trigger a rejection. Instead, we rely solely on the *percentage* filter below.

# --outFilterMismatchNoverReadLmax 0.04
# Limits mismatches to 4% of the length. 
# Example: 100bp read = max 4 mismatches; 200bp pair = max 8 mismatches.

# --outFilterType BySJout
# Only retains junction-crossing reads if the junction is supported by multiple 
# reads or already present in the GTF. Reduces "spurious" (noise) junctions.

# --alignEndsType Local
# Allows STAR to ignore (clip) messy ends or adapter remnants. This "saves" reads 
# from being rejected by the 0.04 mismatch filter if the errors are at the ends.

# --alignIntronMin 20
# Minimum gap to be considered a biological intron. Gaps < 20bp are 
# classified as Deletions (sequencing/mutation errors).

# --alignIntronMax 1000000
# Maximum allowed gap for a single splice (1 million bp). This prevents 
# chimeric-like "jumps" across very large genomic distances.

# --alignMatesGapMax 1000000
# The maximum allowed distance between R1 and R2 mates. Usually set to 
# match the maximum intron size.

# --outSAMunmapped Within
# Includes unmapped reads directly in the main SAM/BAM file with a 
# "0" bitwise flag. This makes it easier to troubleshoot why reads 
# failed to map without needing separate FASTQ files.

# --outSAMtype BAM SortedByCoordinate
# Tells STAR to output a compressed Binary (BAM) file instead of text (SAM) 
# and sorts it by genomic position. This is the required input for 
# indexing (sambamba/samtools) and visualization (IGV).
*/