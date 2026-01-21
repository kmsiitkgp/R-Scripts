// =========================================================================================
// PROCESS: MULTIQC
// =========================================================================================
// Purpose: Aggregates QC reports from multiple tools into a single interactive HTML report
//
// MultiQC scans for output files from supported tools and automatically generates:
//   - Summary statistics tables
//   - Interactive plots
//   - Sample comparison views
//   - Export options (PNG, SVG, CSV)
//
// Supported tools in this pipeline:
//   - FastQC: Base quality, adapter content, duplication
//   - STAR: Alignment rates, multi-mapping statistics
//   - SALMON: Mapping rates, library type detection
//   - RSeQC: Read distribution, gene body coverage, junction stats
//
// The final report allows you to:
//   - Identify failed samples at a glance
//   - Compare QC metrics across samples
//   - Detect batch effects
//   - Make informed decisions about sample exclusion
// =========================================================================================

process MULTIQC {
	
	tag "Merging all QC reports"  // Display in workflow logs
    label 'process_low'  // MultiQC is lightweight (mostly file parsing)
	
	// =================================================================================
	// INPUT
	// =================================================================================
	// Receives a collected list of ALL QC output files from previous processes
	// The channel is created in main.nf by mixing outputs from:
	//   - FASTQC (ZIP files)
	//   - SALMON (quant directories)
	//   - STAR (Log.final.out, ReadsPerGene, SJ.out.tab)
	//   - RSeQC (txt, log, r files)
	input:
	path(all_reports)  // List of all QC files from upstream processes
	
	// =================================================================================
	// OUTPUT
	// =================================================================================
	output:
	path("${params.multiqc_filename()}.html"),  emit: multiqc_html  // Interactive HTML report
	path("${params.multiqc_filename()}_data"),  emit: multiqc_dir   // Data directory (plots, tables)
    path("MULTIQC.error.log"),					emit: multiqc_log   // Process log
	
	// =================================================================================
	// SCRIPT SETUP
	// =================================================================================
	script:	
	
	def LOG = "MULTIQC.error.log"

	"""
	# =============================================================================
	# Run MultiQC to Consolidate All QC Reports
	# =============================================================================
	# MultiQC recursively searches the current directory for recognized output files
	# and aggregates them into a single HTML report
	#
	# Parameters:
	#   --force: Overwrite existing reports
	#   --clean-up: Remove intermediate files after aggregation
	#   --quiet: Suppress progress messages (only show warnings/errors)
	#   --title: Report title (displayed at top of HTML)
	#   --filename: Output filename (without .html extension)
	#   . (dot): Search current directory recursively
	#
	# How MultiQC finds files:
	#   - Nextflow stages all input files in the working directory
	#   - MultiQC has built-in regex patterns for each tool
	#   - It automatically detects file types and parses them
	#   - No need to explicitly list file paths
	#
	# Why use "." instead of "\${all_reports}"?
	#   - Prevents ARG_MAX errors (command line too long)
	#   - Handles directories (like SALMON output) gracefully
	#   - MultiQC is designed to recursively search directories
	#   - More robust to special characters in filenames
	# =============================================================================
	
	multiqc \
		--force \
		--clean-up \
		--quiet \
		--title "${params.multiqc_titlename()}" \
		--filename "${params.multiqc_filename()}" \
		. \
		1>> "${LOG}" 2>&1 \
		|| { echo "❌ ERROR: MultiQC failed" | tee -a "${LOG}" >&2; exit 1; }
	
	echo "✅ SUCCESS: MultiQC completed" >> "${LOG}"		
	
	"""
}

// =========================================================================================
// TECHNICAL NOTES: INPUT HANDLING STRATEGIES
// =========================================================================================
//
// OPTION 1: Use "\${all_reports}" (Explicit File List)
// -----------------------------------------------------
// Pros:
//   - Explicit: You know exactly which files are being processed
//   - Predictable: No surprises from unexpected files in directory
//
// Cons:
//   - ARG_MAX limit: If you have too many files, the command line becomes so long
//     that it exceeds the Linux ARG_MAX (maximum character limit for a command)
//     * Typical limit: ~2MB or ~130,000 characters
//     * For 100 samples with 10 QC files each = 1000 files
//     * Average path length 100 chars = 100,000 chars (approaching limit)
//   - Directory handling: If input includes directories (like SALMON output),
//     the shell may misinterpret the list
//   - Special characters: Spaces or special characters in filenames can cause
//     parsing errors
//
// OPTION 2: Use "." (Current Directory) ✅ RECOMMENDED
// -----------------------------------------------------
// Pros:
//   - No ARG_MAX issues: Short command regardless of file count
//   - Handles directories: MultiQC recursively searches subdirectories
//   - Robust parsing: MultiQC's built-in file detection is reliable
//   - Nextflow staging: All files are already staged in working directory
//
// Cons:
//   - Less explicit: Might process unexpected files if directory is contaminated
//     (unlikely in Nextflow's isolated work directories)
//
// How Nextflow Staging Works:
// ----------------------------
// Before running the script, Nextflow:
//   1. Creates a temporary working directory (work/xx/yyyy...)
//   2. Symlinks every file from input: block into that directory
//   3. Runs the script in that directory
//   4. All outputs are captured from that directory
//
// This means when MultiQC sees ".", it finds all the symlinks Nextflow created
//
// =========================================================================================
// TECHNICAL NOTES: MULTIQC FEATURES
// =========================================================================================
//
// 1. AUTOMATIC TOOL DETECTION:
//    MultiQC recognizes output from 100+ tools using regex patterns:
//    - FastQC: Searches for *_fastqc.zip files
//    - STAR: Searches for Log.final.out files
//    - SALMON: Searches for meta_info.json and quant.sf files
//    - RSeQC: Searches for various .txt, .log files with specific formats
//
// 2. REPORT SECTIONS:
//    The final HTML report contains:
//    a) General Statistics Table:
//       - Sample names
//       - Total reads
//       - Mapping rates
//       - Duplication rates
//       - GC content
//    
//    b) FastQC Section:
//       - Per-base quality plots
//       - Sequence quality distribution
//       - GC content distribution
//       - Adapter content
//       - Duplication levels
//    
//    c) STAR Section:
//       - Uniquely mapped reads
//       - Multi-mapping reads
//       - Unmapped reads
//       - Alignment length distribution
//    
//    d) SALMON Section:
//       - Mapping rates
//       - Library type detection
//       - Fragment length distribution
//    
//    e) RSeQC Section:
//       - Read distribution (CDS/UTR/intron/intergenic)
//       - Gene body coverage
//       - Junction annotation
//       - Insert size distribution
//
// 3. INTERACTIVE FEATURES:
//    - Hover tooltips: See exact values
//    - Sample highlighting: Click to track across plots
//    - Plot export: Download as PNG or SVG
//    - Data export: Download tables as CSV
//    - Plot customization: Toggle samples, change axes
//
// 4. CUSTOMIZATION OPTIONS (add to multiqc command if needed):
//    --config: Custom configuration file
//    --cl-config: Command-line config (e.g., 'sp: {fastqc: {status_checks: false}}')
//    --comment: Add comment to report
//    --sample-names: Rename samples
//    --ignore: Ignore specific file patterns
//    --module: Run specific modules only
//    --exclude: Exclude specific modules
//
// 5. OUTPUT DIRECTORY STRUCTURE:
//    ${params.multiqc_filename()}_data/
//    ├── multiqc_data.json          # All data in JSON format
//    ├── multiqc_general_stats.txt  # General statistics table
//    ├── multiqc_fastqc.txt         # FastQC data
//    ├── multiqc_star.txt           # STAR alignment data
//    ├── multiqc_salmon.txt         # SALMON quantification data
//    ├── multiqc_rseqc_*.txt        # RSeQC data files
//    └── multiqc_sources.txt        # List of source files
//
// 6. COMMON ISSUES AND SOLUTIONS:
//    a) "No analysis results found":
//       - Check that input files are in correct format
//       - Verify file permissions
//       - Check MultiQC log for parsing errors
//    
//    b) Missing samples in report:
//       - Sample names with special characters may cause issues
//       - Check that output files follow expected naming conventions
//       - Look for errors in individual tool outputs
//    
//    c) Very large report file (>50MB):
//       - Many samples (>100) generate large reports
//       - Consider using --data-format to reduce size
//       - Use --flat for simpler plots
//    
//    d) Plots not rendering:
//       - Browser compatibility issues (use Chrome/Firefox)
//       - JavaScript disabled in browser
//       - File corruption during transfer
//
// 7. INTERPRETING THE REPORT:
//    Red flags to look for:
//    - Low mapping rate (<70%): Wrong reference, contamination, or poor quality
//    - High duplication (>50%): Low library complexity or over-amplification
//    - GC bias: PCR amplification artifacts
//    - 3' bias (ratio >3): RNA degradation
//    - High intergenic reads (>10%): Genomic DNA contamination
//    - Adapter content (>5%): Incomplete adapter trimming
//
// 8. NEXT STEPS AFTER MULTIQC:
//    Based on the report, you might:
//    - Exclude failed samples from downstream analysis
//    - Adjust quality filtering thresholds
//    - Perform adapter trimming if not done
//    - Re-extract RNA if quality is poor
//    - Check for batch effects and plan corrections
//    - Document QC decisions for publication
//
// 9. PUBLISHDIR CONFIGURATION:
//    From nextflow.config:
//    - HTML report: ${proj_dir}/06.MultiQC/${project}_MultiQC_Report.html
//    - Data directory: ${proj_dir}/06.MultiQC/${project}_MultiQC_Report_data/
//    - Log file: ${proj_dir}/07.Logs/MULTIQC.error.log
// ========================================================================================='