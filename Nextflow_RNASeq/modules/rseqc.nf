process RSEQC {
	
	def LOG = "${sample_id}.RSEQC.error.log"
	
	input:	
    tuple val(sample_id), path(bam), path(bai)
	val(mode)
	
	output:
	path("${sample_id}*.{pdf,jpeg,png,tiff}"), 	emit: rseqc_plots,		optional: true
    path("${sample_id}*.{txt,log,r}"), 			emit: rseqc_logs, 		optional: true
    path("${LOG}"), 							emit: rseqc_error_log

    script:
	
    // Translate Nextflow 'mode' to RSeQC 'sequencing' flag
    def SEQUENCING_MODE = (mode == "PAIRED_END") ? "PE" : "SE"	
	
	"""
	
	READ_LENGTH=\$(samtools view "${bam}" | head -n1 | awk '{print length(\$10)}')
	
	# 1. Calculate read distribution
	read_distribution.py \
		--input-file "${bam}" \
		--refgene "${params.rseqc_bed}" \
		> "${sample_id}.RSeQC.txt" \
		1>> "${LOG}" 2>&1 \
		|| { echo "❌ ERROR: Read distribution calculation failed for ${sample_id}" | tee -a "${LOG}" >&2; exit 1; }
	
	echo "✅ SUCCESS: Read distribution calculation completed for ${sample_id}" >> "${LOG}"

	
	# 2. Calculate inner distance between read pairs
	if [[ "${mode}" == "PAIRED_END" ]]; then
		inner_distance.py \
			--input-file "${bam}" \
			--refgene "${params.rseqc_bed}" \
			--mapq 30 \
			--out-prefix "${sample_id}" \
			1>> "${LOG}" 2>&1 \
			|| { echo "❌ ERROR: Inner distance calculation failed for ${sample_id}" | tee -a "${LOG}" >&2; exit 1; }
		
		echo "✅ SUCCESS: Inner distance calculation completed for ${sample_id}" >> "${LOG}"

	fi
	
	# 3. Annotate splicing junctions (compared to reference)
	junction_annotation.py \
		--input-file "${bam}" \
		--refgene "${params.rseqc_bed}" \
		--min-intron 50 \
		--mapq 30 \
		--out-prefix "${sample_id}" \
		1>> "${LOG}" 2>&1 \
		|| { echo "❌ ERROR: Junction annotation failed for ${sample_id}" | tee -a "${LOG}" >&2; exit 1; }
	
	echo "✅ SUCCESS: Junction annotation completed for ${sample_id}" >> "${LOG}"		
	
	# 4. Check junction detection saturation (detectability vs sequencing depth)
	junction_saturation.py \
		--input-file "${bam}" \
		--refgene "${params.rseqc_bed}" \
		--min-intron 50 \
		--mapq 30 \
		--out-prefix "${sample_id}" \
		1>> "${LOG}" 2>&1 \
		|| { echo "❌ ERROR: Junction saturation calculation failed for ${sample_id}" | tee -a "${LOG}" >&2; exit 1; }
	
	echo "✅ SUCCESS: Junction saturationn calculation completed for ${sample_id}" >> "${LOG}"

		
	# 5. Profile deletion rates by read position	
	deletion_profile.py \
		--input "${bam}" \
		--read-align-length \${READ_LENGTH} \
		--mapq 30 \
		--out-prefix "${sample_id}" \
		1>> "${LOG}" 2>&1 \
		|| { echo "❌ ERROR: Deletion profile calculation failed for ${sample_id}" | tee -a "${LOG}" >&2; exit 1; }
	
	echo "✅ SUCCESS: Deletion profile calculation completed for ${sample_id}" >> "${LOG}"
		
	
	# 6. Profile mismatch rates by read position
	mismatch_profile.py \
		--input "${bam}" \
		--read-align-length \${READ_LENGTH} \
		--mapq 30 \
		--out-prefix "${sample_id}" \
		1>> "${LOG}" 2>&1 \
		|| { echo "❌ ERROR: Mismatch profile calculation failed for ${sample_id}" | tee -a "${LOG}" >&2; exit 1; }
	
	echo "✅ SUCCESS: Mismatch profile calculation completed for ${sample_id}" >> "${LOG}"
	
	
	# 7. Profile insertion rates by read position	
	insertion_profile.py \
		--input "${bam}" \
		--sequencing ${SEQUENCING_MODE} \
		--mapq 30 \
		--out-prefix "${sample_id}" \
		1>> "${LOG}" 2>&1 \
		|| { echo "❌ ERROR: Insertion profile calculation failed for ${sample_id}" | tee -a "${LOG}" >&2; exit 1; }
	
	echo "✅ SUCCESS: Insertion profile calculation completed for ${sample_id}" >> "${LOG}"
		
	
	# 8. Profile soft-clipping rates at read ends
	clipping_profile.py \
		--input "${bam}" \
		--sequencing ${SEQUENCING_MODE} \
		--mapq 30 \
		--out-prefix "${sample_id}" \
		1>> "${LOG}" 2>&1 \
		|| { echo "❌ ERROR: Clipping profile calculation failed for ${sample_id}" | tee -a "${LOG}" >&2; exit 1; }
	
	echo "✅ SUCCESS: Clipping profile calculation completed for ${sample_id}" >> "${LOG}"
		
	# 9. Calculate coverage uniformity across the gene body (5' to 3' bias)
	geneBody_coverage.py \
		--input "${bam}" \
		--refgene "${params.rseqc_bed}" \
		--minimum_length 100 \
		--format pdf \
		--out-prefix "${sample_id}" \
		1>> "${LOG}" 2>&1 \
		|| { echo "❌ ERROR: Gene body coverage calculation failed for ${sample_id}" | tee -a "${LOG}" >&2; exit 1; }
	
	echo "✅ SUCCESS: Gene body coverage calculation completed for ${sample_id}" >> "${LOG}"	

	"""
}