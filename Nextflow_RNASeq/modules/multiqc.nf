process MULTIQC {
	
	def LOG = "MULTIQC.error.log"
	
	input:
	path(all_reports)

	output:
	path("${params.multiqc_filename}.html"), 	emit: multiqc_html
	path("${params.multiqc_filename}_data"), 	emit: multiqc_dir	
    path("${LOG}"), 							emit: multiqc_log
	
	script:	

	"""
	
	# Consolidate all reports using MultiQC
	multiqc \
		--force \
		--clean-up \
		--quiet \
		--title "${params.multiqc_titlename}" \
		--filename "${params.multiqc_filename}" \
		. \
		1>> "${LOG}" 2>&1 \
		|| { echo "❌ ERROR: MultiQC failed" | tee -a "${LOG}" >&2; exit 1; }
	
	echo "✅ SUCCESS: MultiQC completed" >> "${LOG}"		
	
	"""
}
	
/*
# ==============================================================================================
# You can use ${all_reports} \\ instead of  . \\ but if there are too
# many files, the command becomes so long that it exceeds the Linux ARG_MAX 
# (the maximum character limit for a single command).

# Nextflow is a "staged" environment. Before the script runs, Nextflow creates a
# temporary folder and symlinks every file you listed in your input: block into that folder.
# When you use .: MultiQC scans the current working directory. It sees all the symlinks 
# Nextflow created and processes them.
# When you use "${all_reports}": You are explicitly telling MultiQC exactly which files to
# look at. If for some reason a file name has a space or a special character, the shell 
# might misinterpret the list.

# Remember that your Salmon output is a directory, not a file. If you pass a directory path 
# explicitly in a list, some versions of tools get confused about whether to look at the folder 
# or inside it. MultiQC is designed to "recurse." When it sees . it says "I will search every 
# file and every sub-folder in this directory." This ensures it finds the quant.sf and 
# meta_info.json deep inside your Salmon folders.
# ==============================================================================================
*/
