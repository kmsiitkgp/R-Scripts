workflow VALIDATE_INPUT {

	tag "Checking parameters and file paths"
    
	main:
	
	// =========================================================================================
	// 1. VALIDATE FASTQ NAMING
	// =========================================================================================
	
	// Load all FASTQ files from the raw directory
	def fastq_files_set = files("${params.raw_fastq_dir}/*.f*q.gz")
	def fastq_files_list = fastq_files_set.collect()
	def fastq_files = fastq_files_list.sort()

	// Define validation regex
	def VALID_PATTERN = ~/.*((_Tumor|_Normal))?.*((_R|_r)[12]).*\.f(q|astq)\.gz/
	
	def valid_files   = fastq_files.findAll { it.name ==~ VALID_PATTERN }
	def invalid_files = fastq_files.findAll { !(it.name ==~ VALID_PATTERN) }
		
	// =========================================================================================
	// 2. DETECT FILE FORMAT (fastq.gz vs fq.gz)
	// =========================================================================================
	
	def fq_gz_files     = valid_files.findAll { it.name.endsWith(".fq.gz") }
	def fastq_gz_files  = valid_files.findAll { it.name.endsWith(".fastq.gz") }

	def FILE_FORMAT = ""
	if (fq_gz_files.size() == valid_files.size()) {
		FILE_FORMAT = "fq.gz"
	} else if (fastq_gz_files.size() == valid_files.size()) {
		FILE_FORMAT = "fastq.gz"
	}

	// =========================================================================================
	// 3. DETECT MODE (Strict Case-Sensitive)
	// =========================================================================================
	
	def MODE = ""
	def r1_files = valid_files.findAll { it.name.contains("_r1") }
	def r2_files = valid_files.findAll { it.name.contains("_r2") }
	def R1_files = valid_files.findAll { it.name.contains("_R1") }
	def R2_files = valid_files.findAll { it.name.contains("_R2") }

	if (r1_files.size() == valid_files.size() ||  R1_files.size() == valid_files.size() ) {
		MODE = "SINGLE_END"
	} else if (r1_files.size() > 0 && r1_files.size() == r2_files.size() && r1_files.size() * 2 == valid_files.size()) {
		MODE = "PAIRED_END"
	} else if (R1_files.size() > 0 && R1_files.size() == R2_files.size() && R1_files.size() * 2 == valid_files.size()) {
		MODE = "PAIRED_END"
	}
	
	// =========================================================================================
	// 4. DETECT READ TAGS (For use in filenames/channels)
	// =========================================================================================
	
	def Read1_TAG = ""
	def Read2_TAG = ""

	// Determine R1 Tag
	if (MODE == "SINGLE_END" && r1_files.size() == valid_files.size()) {
		Read1_TAG = "_r1"
	} else if (MODE == "SINGLE_END" && R1_files.size() == valid_files.size()) {
		Read1_TAG = "_R1"
	} else if (MODE == "PAIRED_END" && r1_files.size() * 2 == valid_files.size()){
		Read1_TAG = "_r1"
	} else if (MODE == "PAIRED_END" && R1_files.size() * 2 == valid_files.size()){
		Read1_TAG = "_R1"
	}

	// Determine R2 Tag (Only relevant if MODE is PAIRED_END)
	if (MODE == "PAIRED_END") {
		if (r2_files.size() * 2 == valid_files.size()) {
			Read2_TAG = "_r2"
		} else if (R2_files.size() * 2 == valid_files.size()) {
			Read2_TAG = "_R2"
		}
	}

	// =========================================================================================
	// 5. IMMEDIATE SANITY CHECKS (Fail Fast)
	// =========================================================================================
	
	// A. Check if anything was found at all
	if (valid_files.size() == 0) {
		println "\n" + "!" * 50
		println " ERROR: NO VALID FASTQ FILES FOUND"
		println "!" * 50
		println " Search Path: ${params.raw_fastq_dir}"
		println "!" * 50 + "\n"
		error "Aborting: Zero files matched the pattern in the directory."
	}
	
	// B. Check if files have the wrong names
	if (invalid_files.size() > 0) {		
		println "\n" + "!" * 50
		println " ERROR: ${invalid_files.size()} INVALID FILE(S) DETECTED"
		println "!" * 50		
    
		// Iterates through the list and print just the filenames
		invalid_files.each { 
			println "  - ${it.name}" 
		}    
		println "!" * 50 + "\n"
		error "Aborting: Please rename the files listed above to match VALID_PATTERN."
	} 
		

	// C. Check if all files end with fq.gz or fastq.gz
	if (FILE_FORMAT == "") {
		println "\n" + "!" * 50
		println " FORMAT ERROR: Mixed FASTQ extensions detected"
		println " .fq.gz files    : ${fq_gz_files.size()}"
		println " .fastq.gz files : ${fastq_gz_files.size()}"
		println "!" * 50 + "\n"
		error "Aborting: Please standardize all files to either .fq.gz or .fastq.gz"
	}
	
	// D. Check for Mixed-Case tags (_r1, _r2 OR _R1,_R2)
	if (MODE == ""){
		println "\n" + "!" * 50
		println " ERROR: NAMING VIOLATION"
		println "!" * 50
		println " Lowercase -> r1: ${r1_files.size()}, r2: ${r2_files.size()}"
		println " Uppercase -> R1: ${R1_files.size()}, R2: ${R2_files.size()}"
		println " Total Files Found: ${valid_files.size()}"
		println "!" * 50 + "\n"
		error "Mixed case (_r1 vs _R1) or mismatched pairs detected. Please standardize naming."
	}	
	
	// =========================================================================================
	// 6. SAMPLE COUNTS (Tumor vs Normal)
	// =========================================================================================
	
	def tumor_files  = valid_files.findAll { it.name.contains("_Tumor") }
	def normal_files = valid_files.findAll { it.name.contains("_Normal") }
	
	def TOTAL_SAMPLES    = 0
	def N_TUMOR_SAMPLES  = 0
	def N_NORMAL_SAMPLES = 0	

	if (MODE == "SINGLE_END") {
		TOTAL_SAMPLES    = valid_files.size()
		N_TUMOR_SAMPLES  = tumor_files.size()
		N_NORMAL_SAMPLES = normal_files.size()
	} else if (MODE == "PAIRED_END") {
		// .intdiv(2) ensures we get a whole number (e.g., 10 / 2 = 5)
		TOTAL_SAMPLES    = valid_files.size().intdiv(2)
		N_TUMOR_SAMPLES  = tumor_files.size().intdiv(2)
		N_NORMAL_SAMPLES = normal_files.size().intdiv(2)
	}
	
	// =========================================================================================
	// 7. CREATE SAMPLES CHANNEL
	// =========================================================================================
	
	// Combine all R1 files (handles both r1 and R1 cases)
	def all_r1_files = (r1_files + R1_files).sort()
	def R1_FASTQS_ch = Channel.fromList(all_r1_files)

	def grouped_samples_ch = R1_FASTQS_ch.map { r1 ->
		
		// 1. Extract Sample ID by removing the Read1_TAG
		// If idx is -1, use simpleName (filename minus all extensions) as a fallback		
		def idx = r1.name.lastIndexOf(Read1_TAG)		
		def sample_id = (idx != -1) ? r1.name.take(idx) : r1.simpleName
		
		
		// 2. Pair the files based on MODE
		if (MODE == "PAIRED_END") {
			
			// Use the reverse trick ONLY on the filename to avoid folder-name accidents
			def r1_name = r1.name
			def r2_name = r1_name.reverse().replaceFirst(Read1_TAG.reverse(), Read2_TAG.reverse()).reverse()
			
			// Get the parent directory of r1 and look for r2 there
			def r2 = r1.parent.resolve(r2_name)
			
			// checkIfExists: true is still a great idea!
			if (!r2.exists()) { 
				println "\n" + "!" * 50
				println " ERROR: R2 PAIR NOT FOUND"
				println "!" * 50
				println " Expected R2 : ${r2_name}"
				println " At Location : ${r1.parent}"
				println "!" * 50 + "\n"
				error "Aborting: Missing paired-end mate for ${r1_name}"
			}
					
			return [ sample_id, [r1, r2] ]
			
		} else {
			// SINGLE_END: Return as a list so join(' ') works in your script block
			return [ sample_id, [r1] ]  
		}
	}
	
	// =========================================================================================
	// 8. PRINT SUMMARY STATISTICS
	// =========================================================================================
	
	println "\n" + "=" * 60
	println "              FASTQ INPUT STATISTICS"
	println "=" * 60
	
	// 1. File Format Info
	println " FILE FORMAT      : $FILE_FORMAT"
	println "   - .fq.gz       : ${fq_gz_files.size()}"
	println "   - .fastq.gz    : ${fastq_gz_files.size()}"
	println "-" * 60
	
	// 2. Tag & Mode Info
	println " PIPELINE MODE    : $MODE"
	println " R1 TAG DETECTED  : $Read1_TAG"
	println " R2 TAG DETECTED  : $Read2_TAG"
	println " TAG COUNTS :"
	println "   - Uppercase (R1/R2): ${R1_files.size()} / ${R2_files.size()}"
	println "   - Lowercase (r1/r2): ${r1_files.size()} / ${r2_files.size()}"
	println "-" * 60
	
	// 3. Sample Breakdown
	println " TOTAL SAMPLES        : $TOTAL_SAMPLES"	
	println "   - TUMOR SAMPLES    : $N_TUMOR_SAMPLES"
	println "   - NORMAL SAMPLES   : $N_NORMAL_SAMPLES"
	println " TOTAL FASTQ FILES    : ${valid_files.size()}"
	println "=" * 60 + "\n"
	
	emit:
		// 1. Export Metadata
		mode             = MODE
		read1_tag        = Read1_TAG
		read2_tag        = Read2_TAG
		file_format      = FILE_FORMAT 
		total_samples    = TOTAL_SAMPLES
		n_tumor_samples  = N_TUMOR_SAMPLES
		n_normal_samples = N_NORMAL_SAMPLES
		
		// 2. Export Raw File Channels (Useful for general QC like MultiQC)
		fastq_ch   = Channel.fromList(valid_files)		
		tumor_ch   = Channel.fromList(tumor_files)
		normal_ch  = Channel.fromList(normal_files)

		// 3. Export the "Workhorse" Channel (The grouped tuples)
		samples = grouped_samples_ch
}
