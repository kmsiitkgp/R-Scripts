process SALMON_QUANT {
	
	def LOG = "${sample_id}.SALMON_QUANT.error.log"
	
	input:
	tuple val(sample_id), path(fastq_files)
	path(salmon_index_dir)
	
	output:
    tuple val(sample_id), path("${sample_id}"), 	emit: salmon_quant		
    path("${LOG}"), 								emit: salmon_error_log
	
	script:	
	
	// Determine if paired-end or single-end (Groovy, outside bash)
	def MATES_ARGS = fastq_files.size() == 2 ?
		"--mates1 ${fastq_files[0]} --mates2 ${fastq_files[1]}" :
		"--unmatedReads ${fastq_files[0]}"	
	
	"""
	
	# Salmon quantification
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