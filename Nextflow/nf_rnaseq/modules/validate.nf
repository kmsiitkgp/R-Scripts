// =========================================================================================
// MODULE: VALIDATE_INPUT
// =========================================================================================
// Purpose: Validates FASTQ file naming conventions, detects sequencing mode (SE/PE),
//          and creates organized sample channels for downstream processes.
//
// Key Functions:
//   1. Validates FASTQ naming against expected patterns
//   2. Detects single-end vs paired-end sequencing
//   3. Identifies file format (.fq.gz vs .fastq.gz)
//   4. Groups R1/R2 pairs and creates sample channels
//   5. Provides comprehensive statistics about the input data
// =========================================================================================

workflow VALIDATE_INPUT {	
    
	take:
    data_dir  // String path to directory containing raw FASTQ files
	
	main:
	log.info "==> Validating parameters and file paths"
	
	// =========================================================================================
	// SAFETY CHECK: Prevent null directory errors
	// =========================================================================================
	// This check prevents the /root error that occurs when data_dir is null or empty
	if ( !data_dir ) {
		error "VALIDATE_INPUT: data_dir is null or empty. Check your config logic."
	}

	// Print the resolved path for debugging
	println "DEBUG: data_dir = ${data_dir.toString()}"

	
	// =========================================================================================
	// 1. VALIDATE FASTQ NAMING CONVENTION
	// =========================================================================================
	// Expected naming patterns:
	//   - Optional: _Tumor or _Normal suffix
	//   - Required: _R1/_R2 or _r1/_r2 to indicate read pair
	//   - Extension: .fq.gz or .fastq.gz
	// Examples: 
	//   - Sample1_R1.fastq.gz, Sample1_R2.fastq.gz
	//   - Sample2_Tumor_r1.fq.gz, Sample2_Tumor_r2.fq.gz
	// =========================================================================================
	
	// Convert data_dir to string and glob all compressed FASTQ files
	fastq_dir = data_dir.toString()
    fastq_files_set = files("${fastq_dir}/*.f*q.gz")
	fastq_files_list = fastq_files_set.collect()
	fastq_files = fastq_files_list.sort()

	// Define validation regex pattern
	// Breakdown: .*((_Tumor|_Normal))?.*((_R|_r)[12]).*\.f(q|astq)\.gz
	//   - .* = any characters at start
	//   - ((_Tumor|_Normal))? = optional tumor/normal designation
	//   - .* = any characters in middle
	//   - ((_R|_r)[12]) = required _R1/_R2 or _r1/_r2
	//   - .* = any characters before extension
	//   - \.f(q|astq)\.gz = .fq.gz or .fastq.gz extension
	def VALID_PATTERN = ~/.*((_Tumor|_Normal))?.*((_R|_r)[12]).*\.f(q|astq)\.gz/
	
	// Separate files into valid and invalid based on pattern matching
	valid_files   = fastq_files.findAll { it.name ==~ VALID_PATTERN }
	invalid_files = fastq_files.findAll { !(it.name ==~ VALID_PATTERN) }
		
	// =========================================================================================
	// 2. DETECT FILE FORMAT (.fq.gz vs .fastq.gz)
	// =========================================================================================
	// The pipeline requires all files to use the same extension for consistency
	// Mixed extensions can cause issues in downstream processing
	// =========================================================================================
	
	fq_gz_files     = valid_files.findAll { it.name.endsWith(".fq.gz") }
	fastq_gz_files  = valid_files.findAll { it.name.endsWith(".fastq.gz") }

	def FILE_FORMAT = ""
	if (fq_gz_files.size() == valid_files.size()) {
		FILE_FORMAT = "fq.gz"
	} else if (fastq_gz_files.size() == valid_files.size()) {
		FILE_FORMAT = "fastq.gz"
	}
	// If FILE_FORMAT remains empty, we have mixed formats (error detected later)

	// =========================================================================================
	// 3. DETECT SEQUENCING MODE (SINGLE_END vs PAIRED_END)
	// =========================================================================================
	// Strategy: Check for consistent use of either _r1/_r2 OR _R1/_R2 (case-sensitive)
	// This strict checking prevents mixed-case naming which causes pairing errors
	// =========================================================================================
	
	def MODE = ""
	// Count files with lowercase tags
	r1_files = valid_files.findAll { it.name.contains("_r1") }
	r2_files = valid_files.findAll { it.name.contains("_r2") }
	// Count files with uppercase tags
	R1_files = valid_files.findAll { it.name.contains("_R1") }
	R2_files = valid_files.findAll { it.name.contains("_R2") }

	// Single-end detection: All files are R1 (or r1), no R2 files exist
	if (r1_files.size() == valid_files.size() ||  R1_files.size() == valid_files.size() ) {
		MODE = "SINGLE_END"
	} 
	// Paired-end detection (lowercase): Equal R1/R2 files, and together they account for all files
	else if (r1_files.size() > 0 && r1_files.size() == r2_files.size() && r1_files.size() * 2 == valid_files.size()) {
		MODE = "PAIRED_END"
	} 
	// Paired-end detection (uppercase): Equal R1/R2 files, and together they account for all files
	else if (R1_files.size() > 0 && R1_files.size() == R2_files.size() && R1_files.size() * 2 == valid_files.size()) {
		MODE = "PAIRED_END"
	}
	// If MODE remains empty, we have mixed or mismatched files (error detected later)
	
	// =========================================================================================
	// 4. DETECT READ TAGS (For use in sample ID extraction and pairing)
	// =========================================================================================
	// These tags will be used to:
	//   - Extract sample IDs by removing the tag from filenames
	//   - Pair R1 and R2 files in paired-end mode
	// =========================================================================================
	
	def Read1_TAG = ""
	def Read2_TAG = ""

	// Determine R1 Tag based on detected MODE
	if (MODE == "SINGLE_END" && r1_files.size() == valid_files.size()) {
		Read1_TAG = "_r1"
	} else if (MODE == "SINGLE_END" && R1_files.size() == valid_files.size()) {
		Read1_TAG = "_R1"
	} else if (MODE == "PAIRED_END" && r1_files.size() * 2 == valid_files.size()){
		Read1_TAG = "_r1"
	} else if (MODE == "PAIRED_END" && R1_files.size() * 2 == valid_files.size()){
		Read1_TAG = "_R1"
	}

	// Determine R2 Tag (Only relevant for PAIRED_END mode)
	if (MODE == "PAIRED_END") {
		if (r2_files.size() * 2 == valid_files.size()) {
			Read2_TAG = "_r2"
		} else if (R2_files.size() * 2 == valid_files.size()) {
			Read2_TAG = "_R2"
		}
	}

	// =========================================================================================
	// 5. IMMEDIATE SANITY CHECKS (Fail Fast Approach)
	// =========================================================================================
	// These checks ensure we catch problems early before wasting compute resources
	// Each check provides detailed error messages to help users fix their data
	// =========================================================================================
	
	// CHECK A: Verify at least one valid file was found
	if (valid_files.size() == 0) {
		println "\n" + "!" * 50
		println " ERROR: NO VALID FASTQ FILES FOUND"
		println "!" * 50
		println " Search Path: ${params.raw_fastq_dir()}"
		println "!" * 50 + "\n"
		error "Aborting: Zero files matched the pattern in the directory."
	}
	
	// CHECK B: Report any files with incorrect naming
	if (invalid_files.size() > 0) {		
		println "\n" + "!" * 50
		println " ERROR: ${invalid_files.size()} INVALID FILE(S) DETECTED"
		println "!" * 50		
    
		// List each invalid filename to help user identify issues
		invalid_files.each { 
			println "  - ${it.name}" 
		}    
		println "!" * 50 + "\n"
		error "Aborting: Please rename the files listed above to match VALID_PATTERN."
	} 
		
	// CHECK C: Ensure all files use the same extension
	if (FILE_FORMAT == "") {
		println "\n" + "!" * 50
		println " FORMAT ERROR: Mixed FASTQ extensions detected"
		println " .fq.gz files    : ${fq_gz_files.size()}"
		println " .fastq.gz files : ${fastq_gz_files.size()}"
		println "!" * 50 + "\n"
		error "Aborting: Please standardize all files to either .fq.gz or .fastq.gz"
	}
	
	// CHECK D: Verify consistent case in read pair tags (_r1/_r2 OR _R1/_R2)
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
	// 6. SAMPLE COUNTS (Tumor vs Normal Classification)
	// =========================================================================================
	// Count samples based on _Tumor and _Normal designations in filenames
	// Useful for paired tumor-normal analyses and reporting
	// =========================================================================================
	
	tumor_files  = valid_files.findAll { it.name.contains("_Tumor") }
	normal_files = valid_files.findAll { it.name.contains("_Normal") }
	
	def TOTAL_SAMPLES    = 0
	def N_TUMOR_SAMPLES  = 0
	def N_NORMAL_SAMPLES = 0	

	if (MODE == "SINGLE_END") {
		// In SE mode: 1 file = 1 sample
		TOTAL_SAMPLES    = valid_files.size()
		N_TUMOR_SAMPLES  = tumor_files.size()
		N_NORMAL_SAMPLES = normal_files.size()
	} else if (MODE == "PAIRED_END") {
		// In PE mode: 2 files (R1 + R2) = 1 sample
		// .intdiv(2) ensures integer division (e.g., 10 / 2 = 5)
		TOTAL_SAMPLES    = valid_files.size().intdiv(2)
		N_TUMOR_SAMPLES  = tumor_files.size().intdiv(2)
		N_NORMAL_SAMPLES = normal_files.size().intdiv(2)
	}
	
	// =========================================================================================
	// 7. CREATE SAMPLES CHANNEL
	// =========================================================================================
	// This is the core output: a channel of [sample_id, [fastq_files]] tuples
	// For PE: [sample_id, [R1.fq.gz, R2.fq.gz]]
	// For SE: [sample_id, [R1.fq.gz]]
	// =========================================================================================
	
	// Combine all R1 files (handles both r1 and R1 cases) and sort alphabetically
	all_r1_files = (r1_files + R1_files).sort()
	R1_FASTQS_ch = Channel.fromList(all_r1_files)

	// Map each R1 file to create [sample_id, files] tuples
	grouped_samples_ch = R1_FASTQS_ch.map { r1 ->
		
		// STEP 1: Extract Sample ID by removing the Read1_TAG
		// Example: "Sample1_Tumor_R1.fastq.gz" -> "Sample1_Tumor"
		def idx = r1.name.lastIndexOf(Read1_TAG)
		// If idx is -1 (tag not found), use simpleName as fallback
		def sample_id = (idx != -1) ? r1.name.take(idx) : r1.simpleName
		
		
		// STEP 2: Pair the files based on sequencing MODE
		if (MODE == "PAIRED_END") {
			
			// Use the string reverse trick to safely replace only the R1 tag with R2
			// This avoids accidentally replacing R1 in the sample name
			// Example: "Sample_R1_xyz_R1.fq.gz" -> correctly becomes "Sample_R1_xyz_R2.fq.gz"
			def r1_name = r1.name
			def r2_name = r1_name.reverse().replaceFirst(Read1_TAG.reverse(), Read2_TAG.reverse()).reverse()
			
			// Get the parent directory of r1 and construct the full path to r2
			def r2 = r1.parent.resolve(r2_name)
			
			// Verify that the R2 file actually exists
			if (!r2.exists()) { 
				println "\n" + "!" * 50
				println " ERROR: R2 PAIR NOT FOUND"
				println "!" * 50
				println " Expected R2 : ${r2_name}"
				println " At Location : ${r1.parent}"
				println "!" * 50 + "\n"
				error "Aborting: Missing paired-end mate for ${r1_name}"
			}
			
			// Return tuple: [sample_id, [R1_file, R2_file]]		
			return [ sample_id, [r1, r2] ]
			
		} else {
			// SINGLE_END mode
			// Return as a list so downstream processes can use consistent syntax
			// Return tuple: [sample_id, [R1_file]]
			return [ sample_id, [r1] ]  
		}
	}
	
	// =========================================================================================
	// 8. PRINT SUMMARY STATISTICS
	// =========================================================================================
	// Provides a comprehensive overview of the detected input files
	// Helps users verify the pipeline correctly interpreted their data
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
		// =================================================================================
		// OUTPUT 1: Metadata about the detected files
		// =================================================================================
		// These values provide information about the sequencing run characteristics
		// Downstream processes use these to adapt their behavior (e.g., RSeQC uses mode)
		mode             = MODE             // "SINGLE_END" or "PAIRED_END"
		read1_tag        = Read1_TAG        // "_R1" or "_r1"
		read2_tag        = Read2_TAG        // "_R2" or "_r2" (empty for SE)
		file_format      = FILE_FORMAT      // "fq.gz" or "fastq.gz"
		total_samples    = TOTAL_SAMPLES    // Total number of unique samples
		n_tumor_samples  = N_TUMOR_SAMPLES  // Number of tumor samples
		n_normal_samples = N_NORMAL_SAMPLES // Number of normal samples
		
		// =================================================================================
		// OUTPUT 2: Raw file channels (useful for QC that needs all files)
		// =================================================================================
		fastq_ch   = Channel.fromList(valid_files)    // All valid FASTQ files
		tumor_ch   = Channel.fromList(tumor_files)    // Only tumor FASTQ files
		normal_ch  = Channel.fromList(normal_files)   // Only normal FASTQ files

		// =================================================================================
		// OUTPUT 3: The "workhorse" channel - grouped sample tuples
		// =================================================================================
		// This is the primary output used by alignment/quantification processes
		// Format: [sample_id, [fastq_file(s)]]
		//   - PE example: ["Sample1", [R1.fq.gz, R2.fq.gz]]
		//   - SE example: ["Sample2", [R1.fq.gz]]
		samples = grouped_samples_ch
}