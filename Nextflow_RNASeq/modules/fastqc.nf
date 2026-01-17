process FASTQC {
	
	def LOG = "${sample_id}.${read_type}.FASTQC.error.log"
	
	input:
	tuple val(sample_id), path(fastq_files), val(read_type)
	
	output:
	path("*_fastqc.zip"),	emit: fastqc_zip
	path("*_fastqc.html"),	emit: fastqc_html
	path("${LOG}"),			emit: fastqc_error_log
	
	script:
		
	"""
	
	# Run FastQC to assess quality
	fastqc \
		--threads "${task.cpus}" \
		--quiet \
		${fastq_files.join(' ')} \
		1>> "${LOG}" 2>&1 \
		|| { echo "❌ ERROR: FastQC failed for ${sample_id}" | tee -a "${LOG}" >&2; exit 1; }	
	
	echo "✅ SUCCESS: FastQC completed for ${sample_id}" >> "${LOG}"
	
	"""
}