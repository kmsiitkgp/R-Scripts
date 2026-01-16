process RSEQC {
	
	input:
	// 1. Define 'mode' as an input here
    tuple val(sample_id), path(bam), path(bai), val(mode)

    script:
	
    // 2. Now you can translate Nextflow 'mode' to RSeQC 'sequencing' flag
    def SEQUENCING_MODE = (mode == "PAIRED_END") ? "PE" : "SE"	
	
	"""
	
	READ_LENGTH=\$(samtools view "${bam}" | head -n1 | awk '{print length(\$10)}')
	
	# 1. Calculate read distribution
	read_distribution.py \\
		--input-file "${bam}" \\
		--refgene "${params.rseqc_bed}" \\
		> "${sample_id}.RSeQC.txt" \\
		1>> "${sample_id}.RSEQC.error.log" 2>&1 \\
		&& echo "✅ Read distribution calculation completed successfully." \\
		|| { echo "❌ Read distribution calculation failed. Check ${sample_id}.RSEQC.error.log"; exit 1; }
	
	# 2. Calculate inner distance between read pairs
	if [[ "${mode}" == "PAIRED_END" ]]; then
		inner_distance.py \\
			--input-file "${bam}" \\
			--refgene "${params.rseqc_bed}" \\
			--mapq 30 \\
			--out-prefix "${sample_id}" \\
			1>> "${sample_id}.RSEQC.error.log" 2>&1 \\
			&& echo "✅ Inner distance calculation completed successfully." \\
			|| { echo "❌ Inner distance calculation failed. Check ${sample_id}.RSEQC.error.log"; exit 1; }
	fi
	
	# 3. Annotate splicing junctions (compared to reference)
	junction_annotation.py \\
		--input-file "${bam}" \\
		--refgene "${params.rseqc_bed}" \\
		--min-intron 50 \\
		--mapq 30 \\
		--out-prefix "${sample_id}" \\
		1>> "${sample_id}.RSEQC.error.log" 2>&1 \\
		&& echo "✅ Junction annotation calculation completed successfully." \\
		|| { echo "❌ Junction annotation calculation failed. Check ${sample_id}.RSEQC.error.log"; exit 1; }
	
	# 4. Check junction detection saturation (detectability vs sequencing depth)
	junction_saturation.py \\
		--input-file "${bam}" \\
		--refgene "${params.rseqc_bed}" \\
		--min-intron 50 \\
		--mapq 30 \\
		--out-prefix "${sample_id}" \\
		1>> "${sample_id}.RSEQC.error.log" 2>&1 \\
		&& echo "✅ Junction saturation calculation completed successfully." \\
		|| { echo "❌ Junction saturation calculation failed. Check ${sample_id}.RSEQC.error.log"; exit 1; }
		
	# 5. Profile deletion rates by read position	
	deletion_profile.py \\
		--input "${bam}" \\
		--read-align-length \${READ_LENGTH} \\
		--mapq 30 \\
		--out-prefix "${sample_id}" \\
		1>> "${sample_id}.RSEQC.error.log" 2>&1 \\
		&& echo "✅ Deletion profile calculation completed successfully." \\
		|| { echo "❌ Deletion profile calculation failed. Check ${sample_id}.RSEQC.error.log"; exit 1; }
	
	# 6. Profile mismatch rates by read position
	mismatch_profile.py \\
		--input "${bam}" \\
		--read-align-length \${READ_LENGTH} \\
		--mapq 30 \\
		--out-prefix "${sample_id}" \\
		1>> "${sample_id}.RSEQC.error.log" 2>&1 \\
		&& echo "✅ Mismatch profile calculation completed successfully." \\
		|| { echo "❌ Mismatch profile calculation failed. Check ${sample_id}.RSEQC.error.log"; exit 1; }
	
	# 7. Profile insertion rates by read position	
	insertion_profile.py \\
		--input "${bam}" \\
		--sequencing ${SEQUENCING_MODE} \\
		--mapq 30 \\
		--out-prefix "${sample_id}" \\
		1>> "${sample_id}.RSEQC.error.log" 2>&1 \\
		&& echo "✅ Insertion profile calculation completed successfully." \\
		|| { echo "❌ Insertion profile calculation failed. Check ${sample_id}.RSEQC.error.log"; exit 1; }	
	
	# 8. Profile soft-clipping rates at read ends
	clipping_profile.py \\
		--input "${bam}" \\
		--sequencing ${SEQUENCING_MODE} \\
		--mapq 30 \\
		--out-prefix "${sample_id}" \\
		1>> "${sample_id}.RSEQC.error.log" 2>&1 \\
		&& echo "✅ Clipping profile calculation completed successfully." \\
		|| { echo "❌ Clipping profile calculation failed. Check ${sample_id}.RSEQC.error.log"; exit 1; }
		
	# 9. Calculate coverage uniformity across the gene body (5' to 3' bias)
	geneBody_coverage.py \\
		--input "${bam}" \\
		--refgene "${params.rseqc_bed}" \\
		--minimum_length 100 \\
		--format pdf \\
		--out-prefix "${sample_id}" \\
		1>> "${sample_id}.RSEQC.error.log" 2>&1 \\
		&& echo "✅ Gene body coverage calculation completed successfully." \\
		|| { echo "❌ Gene body coverage calculation failed. Check ${sample_id}.RSEQC.error.log"; exit 1; }
			
	"""
	
	output:
	path("${sample_id}*.{pdf,jpeg,png,tiff}"), 	emit: rseqc_plots, optional: true
    path("${sample_id}*.{txt,log,r}"), 			emit: rseqc_logs
    path("${sample_id}.RSEQC.error.log"), 		emit: rseqc_error_log
}