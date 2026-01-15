workflow VALIDATE_INPUT {
    
	main:

	// 1. Load all FASTQ files from the raw directory
	// =========================================================================================
	// THE "COLLECT" DISTINCTION: Groovy List vs. Nextflow Channel
	//
	// 1. Groovy .collect() (on a List/Collection):
	//    - Returns: A standard Groovy List immediately.
	//    - Usage: Used to transform every item in a list and return the final list.
	//    - Example:
	//          def my_list = [1, 2, 3].collect { it * 2 }
	//          println my_list
	//          // Output: [2, 4, 6]
	//    - Context: Ideal for pre-run validation because you get the list instantly.
	//
	// 2. Nextflow .collect() (on a Channel):
	//    - Returns: A NEW Channel that emits **one single item (pulse)** containing all items in a list.
	//    - Behavior: It waits for all items flowing through the channel to arrive, then emits the full list.
	//    - Example:
	//          ch = Channel.from('file1', 'file2', 'file3')
	//          collected_ch = ch.collect()
	//          collected_ch.subscribe { files_list ->
	//              println files_list
	//          }
	//          // Output: [file1, file2, file3]
	//          // Notice: The channel emits **one single list** after all items have arrived,
	//          //         not each file individually.
	//
	// WHY THIS MATTERS FOR VALIDATION:
	// - If you want to check the number of files, validate names, or abort early,
	//   you need a **standard Groovy List**. Using Nextflow channel collect() here is "lazy"
	//   and will only emit the list **when the workflow actually runs**, which is too late
	//   for pre-run checks.
	// - Therefore, use Groovy methods like `files().findAll { ... }` or `Channel.fromPath(...).collect().toList()`
	//   carefully depending on whether you need immediate access or a reactive stream.
	//
	// =========================================================================================

	def fastq_files = files("${params.raw_fastq_dir}/*.f*q.gz").collect()

	// 2. Define validation regex
	def VALID_PATTERN = ~/.*((_Tumor|_Normal))?.*((_R|_r)[12]).*\.f(q|astq)\.gz/

	// 3. Split into valid/invalid
	// =========================================================================================
	// FILTERING LOGIC: .filter { }
	// The filter() method iterates through the list and keeps only the elements 
	// where the condition inside the curly braces evaluates to 'true'.
	//
	// Method 1: Explicit Variable Name (Best for clarity)
	//    fastq_files.filter { file -> file.name ==~ VALID_PATTERN }
	//
	// Method 2: Implicit 'it' Variable (Shorter syntax)
	//    fastq_files.filter { it.name ==~ VALID_PATTERN }
	//
	// Method 3: The .matches() Method
	//    fastq_files.filter { it.name.matches(VALID_PATTERN) }
	//    - This is a standard Java/Groovy string method.
	//    - Unlike '==~', which is a Groovy operator, .matches() is a function call.
	//
	// TECHNICAL NOTE ON SYNTAX:
	// - ==~      : Groovy "Match" operator. It's concise and optimized for regex objects.
	// - .matches(): Requires the pattern as an argument. It also performs a "strict" 
	//              match (the entire string must match the pattern).
	// - .name    : Ensures we are checking just "sample_R1.fq.gz" and not the full 
	//              system path (/home/user/data/...).
	// =========================================================================================
	
	def valid_files   = fastq_files.findAll { it.name ==~ VALID_PATTERN }
	def invalid_files = fastq_files.findAll { !(it.name ==~ VALID_PATTERN) }

	// 4. Exit if invalid files exist
	// =========================================================================================
	// NOTE: .each vs. .map
	// 
	// 1. .each (Groovy Method):
	//    - Used for SIDE EFFECTS: Printing to console, logging, or throwing errors.
	//    - It iterates over the collection but returns the ORIGINAL collection unchanged.
	//    - In this block, we use .each because we want to trigger a 'stop' signal (error) 
	//      if bad files exist, not transform the files into something else.
	//
	// 2. .map (Nextflow Operator):
	//    - Used for TRANSFORMATION: Changing one data type into another (e.g., File -> Tuple).
	//    - It creates a NEW Channel emitting the results of the transformation.
	//    - Use .map when you need to "prep" data for a downstream process (like STAR).
	//
	// 3. SYNTAX DIFFERENCE:
	//    - .each is usually used on a LIST (after .collect()).
	//    - .map is used on a CHANNEL (to stream items one-by-one).
	// =========================================================================================

	if (invalid_files.size() > 0) {              					// size of the list = number of invalid files
		println "Error: ${invalid_files.size()} invalid file(s) detected:"
		invalid_files.each { file ->                				// iterate over each invalid file
			println "  - $file"
		}
    error "Aborting due to invalid FASTQ filenames."
	}

	// 5. Classify by Tumor/Normal
	def tumor_files  = valid_files.findAll { it.name.contains("_Tumor") }
	def normal_files = valid_files.findAll { it.name.contains("_Normal") }

	// 6. Classify by Read suffix (R1/R2, case-sensitive)
	def r1_files = valid_files.findAll { it.name.contains("_r1") }
	def r2_files = valid_files.findAll { it.name.contains("_r2") }
	def R1_files = valid_files.findAll { it.name.contains("_R1") }
	def R2_files = valid_files.findAll { it.name.contains("_R2") }

	// 7. Detect FILE_FORMAT
	def fq_gz_files     = valid_files.findAll { it.name.endsWith(".fq.gz") }
	def fastq_gz_files  = valid_files.findAll { it.name.endsWith(".fastq.gz") }

	def FILE_FORMAT = ""
	if (fq_gz_files.size() == valid_files.size()) {
		FILE_FORMAT = "fq.gz"
	} else if (fastq_gz_files.size() == valid_files.size()) {
		FILE_FORMAT = "fastq.gz"
	} else {
		error "Mixed FASTQ formats detected (.fq.gz vs .fastq.gz). Please standardize."
	}

	// 8. Detect MODE
	def MODE = ""
	if (r1_files.size() == valid_files.size() ||  R1_files.size() == valid_files.size() ) {
		MODE = "SINGLE_END"
	} else if (r1_files.size() == r2_files.size() && r1_files.size() == valid_files.size()/2) {
		MODE = "PAIRED_END"
	} else if (R1_files.size() == R2_files.size() && R1_files.size() == valid_files.size()/2) {
		MODE = "PAIRED_END"
	} else {
		error "FASTQ file counts are inconsistent: cannot determine SINGLE_END vs PAIRED_END."
	}

	// Detect READ_SUFFIX
	def READ_SUFFIX = ""
	if (r1_files.size() == valid_files.size() ||  r1_files.size() == valid_files.size()/2 ) {
		READ_SUFFIX = "_r1"
	} else if (R1_files.size() == valid_files.size() || R1_files.size() == valid_files.size()/2) {
		READ_SUFFIX = "_R1"
	} 

	// Detect READ_REPL
	def READ_REPL = ""
	if (r2_files.size() == valid_files.size() ||  r2_files.size() == valid_files.size()/2 ) {
		READ_REPL = "_r2"
	} else if (R2_files.size() == valid_files.size() || R2_files.size() == valid_files.size()/2) {
		READ_REPL = "_R2"
	} 

	// Determine N_TUMOR_SAMPLES and N_NORMAL_SAMPLES
	def N_SAMPLES = 0
	def N_TUMOR_SAMPLES = 0
	def N_NORMAL_SAMPLES = 0
	if (MODE == "SINGLE_END") {
		N_SAMPLES = valid_files.size()
		N_TUMOR_SAMPLES = tumor_files.size()
		N_NORMAL_SAMPLES = normal_files.size()
	} else if (MODE == "PAIRED_END") {
		N_SAMPLES = valid_files.size() / 2
		N_TUMOR_SAMPLES = tumor_files.size() / 2
		N_NORMAL_SAMPLES = normal_files.size() / 2
	} 
	
	// Create the grouped SAMPLES channel
	// We take the R1 files and 'map' them to a tuple: [ SampleID, [R1, R2] ]
	// =========================================================================================
	// SCOPING WITH 'def':
	// We use 'def' to declare variables as 'local' to this specific closure (the code block).
	// This is critical because:
	//   1. It prevents 'leaking' variables into the global pipeline scope.
	//   2. It ensures that parallel tasks don't accidentally overwrite each other's data.
	
	// FILE PATH LOGIC:
	// 1. r1.name     : Returns ONLY the filename (e.g., "sample_R1.fq.gz"). 
	//                  Use this for IDs or renaming.
	// 2. r1.toString(): Returns the ABSOLUTE PATH (e.g., "/work/data/sample_R1.fq.gz").
	//                  REQUIRED when locating the 'partner' R2 file in the same directory.
	
	// THE lastIndexOf() METHOD:
	// We use lastIndexOf instead of indexOf to ensure we find the suffix 
	// closest to the file extension. This prevents issues if the suffix 
	// (e.g., "R1") appears earlier in the sample name itself.

	// THE take(idx) METHOD:
	// This Groovy method slices the string from the beginning (index 0) up to the specified
	// index (idx). It effectively "chops off" the suffix and everything that follows it 
	// (like .fastq.gz), leaving only the prefix as the Sample ID.
	
	// THE file() HELPER:
	// This converts a raw string (the path we just built) into a formal Nextflow 'File Object'.
	// This is critical because:
	//   a) It allows Nextflow to create the symlink 'shortcuts' in the work directory.
	//   b) Using 'checkIfExists: true' ensures the pipeline crashes early if R2 is missing.
	// =========================================================================================
		
	def R1_FASTQS_ch = Channel.fromList((r1_files + R1_files).sort())
	def grouped_samples_ch = R1_FASTQS_ch.map { r1 ->
		// Extract Sample ID
		def idx = r1.name.lastIndexOf(READ_SUFFIX)
		def sample_id = (idx != -1) ? r1.name.take(idx) : r1.simpleName
		
		// Create tuple
		if (MODE == "PAIRED_END") {
			def r2_path = r1.toString().reverse().replaceFirst(READ_SUFFIX.reverse(), READ_REPL.reverse()).reverse()
			def r2 = file(r2_path, checkIfExists: true)
			return [ sample_id, [r1, r2] ]
			
		} else if (MODE == "SINGLE_END"){			
			return [ sample_id, [r1] ]  
			// keep r1 as list so STAR --readFilesIn ${fastq_files.join(' ')} doesnt fail for single end
		}
	}

	// 11. Mixed-case check for read suffix
	if ((R1_files.size()>0 && r1_files.size()>0) || 
		(R2_files.size()>0 && r2_files.size()>0) || 
		(R1_files.size()>0 && r2_files.size()>0) ||
		(R2_files.size()>0 && r1_files.size()>0)) {
			error "Mixed-case detected for read suffixes! Use ALL uppercase (_R1/_R2) or ALL lowercase (_r1/_r2)."
	}

	// 12. Basic sanity checks
	if (valid_files.size()  <= 0) {
		error "No FASTQ files found in ${params.raw_fastq_dir}."
	}

	if (N_SAMPLES <= 0) {
		error "Number of samples is zero. Check your R1/R2 counts."
	}

	if (MODE == "PAIRED_END" && valid_files.size() % 2 != 0) {
		error "Paired-end mode requires an even number of FASTQ files."
	}

	// 13. Print stats
	if (params.print_stats) {
		println "=========================================================="
		println "--- FASTQ Input Statistics ---"    
		
		println "fq.gz Files			: ${fq_gz_files.size()}"
		println "fastq.gz Files			: ${fastq_gz_files.size()}"		
		println "File format            : $FILE_FORMAT" 	
		
		println "R1 Files				: ${R1_files.size()}"
		println "R2 Files				: ${R2_files.size()}"
		println "r1 Files				: ${r1_files.size()}"
		println "r2 Files				: ${r2_files.size()}"
		println "Read suffix            : $READ_SUFFIX"
		println "Read replacement       : $READ_REPL"
		
		println "Mode                   : $MODE"
		println "Total FASTQ files      : ${valid_files.size()}"	
		println "Total samples    	    : $N_SAMPLES"
		
		println "Tumor samples          : $N_TUMOR_SAMPLES"
		println "Normal samples         : $N_NORMAL_SAMPLES"
		
		println "=========================================================="
	}

	emit:
		// Export the variables so main.nf can see them
		mode             = MODE
		read_suffix      = READ_SUFFIX
		read_repl        = READ_REPL
		file_format 	 = FILE_FORMAT 
		n_samples  		 = N_SAMPLES
		n_tumor_samples  = N_TUMOR_SAMPLES
		n_normal_samples = N_NORMAL_SAMPLES
		
		// Export the channels so main.nf can see them
		fastq_ch   = Channel.fromList(valid_files)		
		tumor_ch   = Channel.fromList(tumor_files)
		normal_ch  = Channel.fromList(normal_files)
		samples_ch = grouped_samples_ch 
