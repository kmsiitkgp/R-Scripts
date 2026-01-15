process FASTQC {

	input:
	path(fastq_file)
	
	output:
	
	script:
	
	"""
	
	fastqc \\
		--threads "${task.cpus}" \\
		--quiet \\
		"${fastq_file}" \\
		1>> "${fastq_file}.FASTQC.error.log" 2>&1 \\
		&& echo "✅ FastQC completed successfully." \\
		|| { echo "❌ FastQC failed. Check ${fastq_file}.FASTQC.error.log"; exit 1; }	
	
	
	"""
	
	
	



















}