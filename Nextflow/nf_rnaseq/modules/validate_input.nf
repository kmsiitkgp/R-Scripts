// =========================================================================================
// WORKFLOW: VALIDATE_INPUT
// =========================================================================================
// Purpose: Validates FASTQ naming, detects SE/PE mode, creates sample channels
//
// What it does:
//   1. Validates FASTQ file naming conventions
//   2. Detects single-end vs paired-end sequencing
//   3. Ensures consistent file format (.fq.gz vs .fastq.gz)
//   4. Pairs R1/R2 files correctly
//   5. Creates organized sample channels for downstream processes
//
// Why this matters:
//   - Catches naming errors BEFORE wasting compute time
//   - Prevents cryptic pairing errors later in pipeline
//   - Provides clear, actionable error messages
//
// For detailed explanation, see: docs/validate_input.md
// =========================================================================================

workflow VALIDATE_INPUT {

    take:
    data_dir                          // Path to directory containing raw FASTQ files

    main:
    log.info "==> Validating parameters and file paths"

    // Safety check: Prevent null directory errors
    if ( !data_dir ) {
        error "VALIDATE_INPUT: data_dir is null or empty. Check your config logic."
    }

    println "DEBUG: data_dir = ${data_dir.toString()}"

    // =================================================================================
    // 1. FIND AND VALIDATE FASTQ FILES
    // =================================================================================

    // Glob all compressed FASTQ files
    fastq_dir = data_dir.toString()
    fastq_files_set = files("${fastq_dir}/*.f*q.gz")
    fastq_files_list = fastq_files_set.collect()
    fastq_files = fastq_files_list.sort()

    // Define validation regex pattern
    // Expected: *_R1.fq.gz, *_R2.fq.gz (or _r1/_r2)
    // Optional: _Tumor or _Normal designation
    def VALID_PATTERN = ~/.*((_Tumor|_Normal))?.*((_R|_r)[12]).*\.f(q|astq)\.gz/

    // Separate valid from invalid files
    valid_files   = fastq_files.findAll { it.name ==~ VALID_PATTERN }
    invalid_files = fastq_files.findAll { !(it.name ==~ VALID_PATTERN) }

    // =================================================================================
    // 2. DETECT FILE FORMAT
    // =================================================================================

    fq_gz_files     = valid_files.findAll { it.name.endsWith(".fq.gz") }
    fastq_gz_files  = valid_files.findAll { it.name.endsWith(".fastq.gz") }

    def FILE_FORMAT = ""
    if (fq_gz_files.size() == valid_files.size()) {
        FILE_FORMAT = "fq.gz"
    } else if (fastq_gz_files.size() == valid_files.size()) {
        FILE_FORMAT = "fastq.gz"
    }

    // =================================================================================
    // 3. DETECT SEQUENCING MODE (SE vs PE)
    // =================================================================================

    def MODE = ""
    r1_files = valid_files.findAll { it.name.contains("_r1") }
    r2_files = valid_files.findAll { it.name.contains("_r2") }
    R1_files = valid_files.findAll { it.name.contains("_R1") }
    R2_files = valid_files.findAll { it.name.contains("_R2") }

    // SE: All files are R1, no R2
    if (r1_files.size() == valid_files.size() ||  R1_files.size() == valid_files.size() ) {
        MODE = "SINGLE_END"
    }
    // PE (lowercase): Equal R1/R2 files
    else if (r1_files.size() > 0 && r1_files.size() == r2_files.size() && r1_files.size() * 2 == valid_files.size()) {
        MODE = "PAIRED_END"
    }
    // PE (uppercase): Equal R1/R2 files
    else if (R1_files.size() > 0 && R1_files.size() == R2_files.size() && R1_files.size() * 2 == valid_files.size()) {
        MODE = "PAIRED_END"
    }

    // =================================================================================
    // 4. DETECT READ TAGS
    // =================================================================================

    def Read1_TAG = ""
    def Read2_TAG = ""

    if (MODE == "SINGLE_END" && r1_files.size() == valid_files.size()) {
        Read1_TAG = "_r1"
    } else if (MODE == "SINGLE_END" && R1_files.size() == valid_files.size()) {
        Read1_TAG = "_R1"
    } else if (MODE == "PAIRED_END" && r1_files.size() * 2 == valid_files.size()){
        Read1_TAG = "_r1"
    } else if (MODE == "PAIRED_END" && R1_files.size() * 2 == valid_files.size()){
        Read1_TAG = "_R1"
    }

    if (MODE == "PAIRED_END") {
        if (r2_files.size() * 2 == valid_files.size()) {
            Read2_TAG = "_r2"
        } else if (R2_files.size() * 2 == valid_files.size()) {
            Read2_TAG = "_R2"
        }
    }

    // =================================================================================
    // 5. VALIDATION CHECKS (Fail Fast)
    // =================================================================================

    // CHECK A: At least one valid file found
    if (valid_files.size() == 0) {
        println "\n" + "!" * 50
        println " ERROR: NO VALID FASTQ FILES FOUND"
        println "!" * 50
        println " Search Path: ${params.raw_fastq_dir()}"
        println " Pattern: *.f*q.gz"
        println " Files found: ${fastq_files.size()}"
        println "!" * 50
        println "\n Possible reasons:"
        println "   1. Wrong directory path in config"
        println "   2. Files not gzipped (.gz extension missing)"
        println "   3. Files don't match naming convention"
        println "   4. Directory is empty"
        println "!" * 50 + "\n"
        error "Aborting: Zero files matched the pattern in the directory."
    }

    // CHECK B: Report invalid file names
    if (invalid_files.size() > 0) {
        println "\n" + "!" * 50
        println " ERROR: ${invalid_files.size()} INVALID FILE(S) DETECTED"
        println "!" * 50
        println " These files don't match the expected naming pattern:"
        println " Expected: *_R1.fq.gz / *_R2.fq.gz (or _r1/_r2)"
        println " Optional: *_Tumor_* or *_Normal_* designation"
        println ""

        invalid_files.each {
            println "  ✗ ${it.name}"
        }

        println ""
        println " Common fixes:"
        println "   1. Rename _1.fq.gz → _R1.fq.gz (add 'R')"
        println "   2. Add .gz extension if missing"
        println "   3. Ensure consistent case (_R1 not _r1 mixed)"
        println "!" * 50 + "\n"
        error "Aborting: Please rename the files listed above to match expected pattern."
    }

    // CHECK C: Ensure consistent file format
    if (FILE_FORMAT == "") {
        println "\n" + "!" * 50
        println " FORMAT ERROR: Mixed FASTQ extensions detected"
        println "!" * 50
        println " .fq.gz files    : ${fq_gz_files.size()}"
        println " .fastq.gz files : ${fastq_gz_files.size()}"
        println ""
        println " All files must use the SAME extension."
        println " Choose one and rename all files:"
        println "   Option 1: Rename all to .fq.gz"
        println "   Option 2: Rename all to .fastq.gz"
        println "!" * 50 + "\n"
        error "Aborting: Please standardize all files to either .fq.gz or .fastq.gz"
    }

    // CHECK D: Verify consistent read pair tags
    if (MODE == ""){
        println "\n" + "!" * 50
        println " ERROR: INCONSISTENT READ PAIR NAMING"
        println "!" * 50
        println " File counts by tag:"
        println "   Lowercase → r1: ${r1_files.size()}, r2: ${r2_files.size()}"
        println "   Uppercase → R1: ${R1_files.size()}, R2: ${R2_files.size()}"
        println "   Total files   : ${valid_files.size()}"
        println ""
        println " Possible issues:"
        println "   1. Mixed case (_r1 and _R1 in same dataset)"
        println "   2. Mismatched pairs (different number of R1 vs R2)"
        println "   3. Missing R2 files for some samples"
        println ""
        println " Solution:"
        println "   - Use ONLY _R1/_R2 (or ONLY _r1/_r2)"
        println "   - Ensure every R1 has matching R2 with same sample name"
        println "!" * 50 + "\n"
        error "Mixed case or mismatched pairs detected. Please standardize naming."
    }

    // =================================================================================
    // 6. COUNT SAMPLES (Tumor/Normal classification)
    // =================================================================================

    tumor_files  = valid_files.findAll { it.name.contains("_Tumor") }
    normal_files = valid_files.findAll { it.name.contains("_Normal") }

    def TOTAL_SAMPLES    = 0
    def N_TUMOR_SAMPLES  = 0
    def N_NORMAL_SAMPLES = 0

    if (MODE == "SINGLE_END") {
        TOTAL_SAMPLES    = valid_files.size()
        N_TUMOR_SAMPLES  = tumor_files.size()
        N_NORMAL_SAMPLES = normal_files.size()
    } else if (MODE == "PAIRED_END") {
        TOTAL_SAMPLES    = valid_files.size().intdiv(2)
        N_TUMOR_SAMPLES  = tumor_files.size().intdiv(2)
        N_NORMAL_SAMPLES = normal_files.size().intdiv(2)
    }

    // =================================================================================
    // 7. CREATE SAMPLE CHANNELS
    // =================================================================================

    all_r1_files = (r1_files + R1_files).sort()
    R1_FASTQS_ch = Channel.fromList(all_r1_files)

    grouped_samples_ch = R1_FASTQS_ch.map { r1 ->

        // Extract sample ID by removing Read1_TAG
        def idx = r1.name.lastIndexOf(Read1_TAG)
        def sample_id = (idx != -1) ? r1.name.take(idx) : r1.simpleName

        if (MODE == "PAIRED_END") {
            // Find R2 mate using reverse-replace trick
            // Reverses string, replaces first (=last) occurrence, reverses back
            def r1_name = r1.name
            def r2_name = r1_name.reverse().replaceFirst(Read1_TAG.reverse(), Read2_TAG.reverse()).reverse()
            def r2 = r1.parent.resolve(r2_name)

            // Verify R2 exists
            if (!r2.exists()) {
                println "\n" + "!" * 50
                println " ERROR: MISSING R2 PAIR"
                println "!" * 50
                println " R1 file  : ${r1.name}"
                println " Expected : ${r2_name}"
                println " Location : ${r1.parent}"
                println ""
                println " This R1 file has no matching R2 file."
                println " Check for:"
                println "   1. Typo in R2 filename"
                println "   2. Missing file (incomplete download/transfer)"
                println "   3. Case mismatch (R1 vs r1)"
                println "!" * 50 + "\n"
                error "Aborting: Missing paired-end mate for ${r1.name}"
            }

            return [ sample_id, [r1, r2] ]

        } else {
            // Single-end: return R1 only (as list for consistency)
            return [ sample_id, [r1] ]
        }
    }

    // =================================================================================
    // 8. PRINT SUMMARY
    // =================================================================================

    println "\n" + "=" * 60
    println "              FASTQ INPUT VALIDATION SUMMARY"
    println "=" * 60

    println " FILE FORMAT      : $FILE_FORMAT"
    println "   - .fq.gz       : ${fq_gz_files.size()}"
    println "   - .fastq.gz    : ${fastq_gz_files.size()}"
    println "-" * 60

    println " SEQUENCING MODE  : $MODE"
    println " READ TAGS USED   :"
    println "   - R1 tag       : $Read1_TAG"
    if (MODE == "PAIRED_END") {
        println "   - R2 tag       : $Read2_TAG"
    }
    println " TAG DISTRIBUTION :"
    println "   - Uppercase (R1/R2): ${R1_files.size()} / ${R2_files.size()}"
    println "   - Lowercase (r1/r2): ${r1_files.size()} / ${r2_files.size()}"
    println "-" * 60

    println " SAMPLE SUMMARY       :"
    println "   - Total samples    : $TOTAL_SAMPLES"
    println "   - Tumor samples    : $N_TUMOR_SAMPLES"
    println "   - Normal samples   : $N_NORMAL_SAMPLES"
    println "   - Other samples    : ${TOTAL_SAMPLES - N_TUMOR_SAMPLES - N_NORMAL_SAMPLES}"
    println " TOTAL FASTQ FILES    : ${valid_files.size()}"
    println "=" * 60 + "\n"

    emit:
        // Metadata outputs
        mode             = MODE
        read1_tag        = Read1_TAG
        read2_tag        = Read2_TAG
        file_format      = FILE_FORMAT
        total_samples    = TOTAL_SAMPLES
        n_tumor_samples  = N_TUMOR_SAMPLES
        n_normal_samples = N_NORMAL_SAMPLES

        // Raw file channels
        fastq_ch   = Channel.fromList(valid_files)
        tumor_ch   = Channel.fromList(tumor_files)
        normal_ch  = Channel.fromList(normal_files)

        // Grouped sample tuples: [sample_id, [R1] or [R1, R2]]
        samples = grouped_samples_ch
}

// =========================================================================================
// QUICK REFERENCE
// =========================================================================================
//
// Valid naming patterns:
//   Sample1_R1.fastq.gz, Sample1_R2.fastq.gz
//   Sample1_Tumor_r1.fq.gz, Sample1_Tumor_r2.fq.gz
//   Patient001_Normal_R1.fq.gz, Patient001_Normal_R2.fq.gz
//
// Invalid patterns:
//   Sample1_1.fq.gz (missing "R")
//   Sample1.fastq (not gzipped)
//   Sample1_R3.fq.gz (should be R1 or R2)
//
// Output channels:
//   samples: [sample_id, [fastq_files]] - main output for downstream
//   mode: "SINGLE_END" or "PAIRED_END"
//   fastq_ch: All valid FASTQ files
//
// For comprehensive guide, see: docs/validate_input.md
// ========================================================================================='