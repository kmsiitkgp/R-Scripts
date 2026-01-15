process SALMON_QUANT {

	input:
	tuple val(sample_id), path(fastq_files)
	
	script:	
	
	// Determine if paired-end or single-end (Groovy, outside bash)
	def MATES_ARGS = fastq_files.size() == 2 ?
		"--mates1 ${fastq_files[0]} --mates2 ${fastq_files[1]}" :
		"--mates1 ${fastq_files[0]}"
	
	"""
	# 1. Salmon quantification
	salmon quant \\
		--validateMappings \\
		--index "${params.salmon_index_dir}" \\
		${MATES_ARGS} \\
		${params.SALMON_ARGS.join(' ')} \\
		--threads "${task.cpus}" \\
		--output "${sample_id}" \\
		1>> "${sample_id}.SALMON.error.log" 2>&1 \\
		&& echo "✅ Salmon quantification completed successfully." \\
		|| { echo "❌ Salmon quantification failed. Check ${sample_id}.SALMON.error.log"; exit 1; }
	
	
	"""
	
	output:
	// 1. Salmon results directory
    tuple val(sample_id), path("${sample_id}"), emit: salmon_quant
	
    // 2. Your custom error log
    path "${sample_id}.SALMON.error.log", emit: salmon_error_log
}