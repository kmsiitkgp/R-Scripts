process CELLRANGER_COUNT {
	
	tag "Running CellRanger count on ${sample_id}"
	label 'process_high'
	
	input:
	val(sample_id)	
	
	output:
	path("${sample_id}/outs/raw_feature_bc_matrix"), 		emit: raw_matrix_dir
	path("${sample_id}/outs/filtered_feature_bc_matrix"), 	emit: filt_matrix_dir
	path("${sample_id}/outs/analysis"), 					emit: analysis_dir
	path("${sample_id}/outs/*.bam*"), 						emit: bam_files, 			optional: true
	path("${sample_id}/outs/*.csv"), 						emit: metric_summary
	path("${sample_id}/outs/*.html"), 						emit: web_summary
	path("${sample_id}/outs/*.h5"), 						emit: h5_files
	path("${sample_id}/outs/*.json"), 						emit: run_info
	path("${sample_id}/outs/*.cloupe"), 					emit: cloupe
	path("*.CELLRANGER_COUNT.error.log"),					emit: cellranger_count_error_log  // Process log	
	
	script:
	
	def LOG = "${sample_id}.CELLRANGER_COUNT.error.log"
	
	"""
	
	cellranger count \
		${params.CELLRANGER_ARGS.join(' ')} \		
		--transcriptome "${params.cellranger_index_dir}" \
		--fastqs "${params.raw_fastq_dir}" \
		--sample "${sample_id}" \
		--id "${sample_id}" \
		--create-bam "${params.create_bam}" \		
		--localcores "${task.cpus}" \
		1>> "${LOG}" 2>&1 \
        || { echo "❌ ERROR: CellRanger count failed for ${sample_id}" | tee -a "${LOG}" >&2; exit 1; }
	
	echo "✅ SUCCESS: CellRanger count completed for ${sample_id}" >> "${LOG}"	
	"""
}

	