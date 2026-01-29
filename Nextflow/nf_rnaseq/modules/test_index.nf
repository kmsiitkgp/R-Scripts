// =========================================================================================
// PROCESS: TEST_INDEX (DEBUGGING UTILITY)
// =========================================================================================
// Purpose: Debugging tool to diagnose Nextflow resume/cache failures
//
// What it does:
//   - Prints all input parameters that affect process hashing
//   - Helps identify which parameter breaks -resume functionality
//   - Tests argument passing strategies (closures vs pre-joined values)
//
// Historical context:
//   This process was created to debug why SALMON_QUANT and STAR_ALIGN failed to resume
//   while FASTQC resumed properly. Root cause: params.SALMON_ARGS().join(' ') inside
//   the process was regenerating different hashes on each run, breaking cache lookup.
//
// Solution discovered:
//   Join arguments in main.nf ONCE and pass as a val() parameter. This ensures
//   consistent hashing across pipeline runs, enabling proper -resume behavior.
//
// When to use:
//   - Debugging resume failures in new processes
//   - Testing parameter passing strategies
//   - Verifying cache key stability
//
// For Nextflow caching mechanics, see: docs/nextflow_caching.md
// =========================================================================================

process TEST_INDEX {

    tag "Testing ${sample_id}"
    // No label needed - this is a debug process with minimal resources

    // =================================================================================
    // INPUT
    // =================================================================================
    input:
    tuple val(sample_id), path(fastq_files)  // Sample identifier and FASTQ files
    path salmon_index_dir                    // Salmon index directory path
    val s_args                               // Pre-joined arguments from main.nf

    // =================================================================================
    // EXECUTION
    // =================================================================================
    script:

    // Recreate the exact logic used in real processes
    // This helps verify that parameter construction doesn't break caching
    def MATES_ARGS = fastq_files.size() == 2 ?
        "--mates1 ${fastq_files[0]} --mates2 ${fastq_files[1]}" :
        "--unmatedReads ${fastq_files[0]}"

    """
    # Print all parameters that contribute to process hash
    # If any of these change between runs, cache will miss
    echo "=== CACHE DEBUG ==="
    echo "Sample ID       : ${sample_id}"
    echo "Index Path      : ${salmon_index_dir}"
    echo "Mates Arguments : ${MATES_ARGS}"
    echo "CPUs            : ${task.cpus}"
    echo "Joined Args     : ${s_args}"
    echo "==================="

    # Explanation:
    # - sample_id: Should be stable (from file parsing)
    # - salmon_index_dir: Should be stable (file path)
    # - MATES_ARGS: Derived from fastq_files, should be stable
    # - task.cpus: Stable if config unchanged
    # - s_args: CRITICAL - must be pre-joined in main.nf to stay stable
    #
    # If s_args is computed inside process (e.g., params.X().join(' ')),
    # it regenerates on each run with different internal IDs, breaking cache.
    """
}

// =========================================================================================
// DEBUGGING GUIDE
// =========================================================================================
//
// How Nextflow caching works:
//   - Each process execution generates a unique hash from:
//     * Input values (val, path)
//     * Script content
//     * Container/environment
//     * Configuration directives
//   - Hash is used to lookup previous results in work/ directory
//   - If hash matches, Nextflow reuses cached outputs (-resume)
//
// Common cache-breaking mistakes:
//   1. Closures inside process (e.g., params.X().join(' '))
//      → Solution: Pre-compute in main.nf, pass as val()
//
//   2. Random/timestamp in script
//      → Solution: Use deterministic values or external files
//
//   3. Changed config between runs
//      → Solution: Keep process-specific configs stable
//
//   4. Dynamic file paths with timestamps
//      → Solution: Use consistent naming schemes
//
// How to use this debug process:
//   1. Add TEST_INDEX to your workflow
//   2. Run pipeline once: nextflow run main.nf
//   3. Run again with -resume: nextflow run main.nf -resume
//   4. Check if TEST_INDEX shows "Cached" or re-runs
//   5. Compare printed values between runs
//   6. Identify which parameter changed unexpectedly
//
// Example workflow integration:
//   ```groovy
//   TEST_INDEX(
//       samples_ch,              // From VALIDATE_INPUT
//       salmon_index_dir,        // From SALMON_INDEX or config
//       params.SALMON_ARGS().join(' ')  // Pre-join BEFORE passing!
//   )
//   ```
//
// Expected output (if caching works):
//   First run:  "Submitted process > TEST_INDEX (Sample1)"
//   Second run: "Cached process > TEST_INDEX (Sample1)"
//
// If it re-runs instead of caching:
//   - One of the echoed values is changing between runs
//   - Most common culprit: "Joined Args" if computed inside process
//
// For comprehensive Nextflow caching guide, see: docs/nextflow_caching.md
// =========================================================================================
