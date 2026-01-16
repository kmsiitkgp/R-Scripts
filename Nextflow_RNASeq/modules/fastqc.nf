process FASTQC {

	input:
	tuple val(sample_id), path(fastq_files)
		
	script:
	
	"""
	# 1. Run FastQC to assess quality
	fastqc \\
		--threads "${task.cpus}" \\
		--quiet \\
		${fastq_files.join(' ')} \\
		1>> "${sample_id}.FASTQC.error.log" 2>&1 \\
		&& echo "✅ FastQC completed successfully." \\
		|| { echo "❌ FastQC failed. Check ${sample_id}.FASTQC.error.log"; exit 1; }	
	
	
	"""
	
	output:
	// MultiQC Inputs
	path "*fastqc.zip", \
		emit: fastqc_zip
	path "*fastqc.html", \
		emit: fastqc_html
	
	// Troubleshooting
	path "${sample_id}.FASTQC.error.log", \
		emit: fastqc_error_log

}