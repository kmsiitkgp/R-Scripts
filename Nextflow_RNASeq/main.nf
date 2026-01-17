#!/usr/bin/env nextflow
nextflow.enable.dsl=2 

include { VALIDATE_INPUT } 				from './modules/validate.nf'
include { FASTQC as FASTQC_RAW } 		from './modules/fastqc.nf'
include { FASTQC as FASTQC_TRIMMED } 	from './modules/fastqc.nf'
include { PREP_REFERENCE } 				from './modules/prep_reference.nf'
include { SALMON_QUANT } 				from './modules/salmon_quant.nf'
include { STAR_ALIGN } 					from './modules/star_align.nf'
include { RSEQC } 						from './modules/rseqc.nf'
include { MULTIQC } 					from './modules/multiqc.nf'

workflow {
	
	// Validate fastqs
    VALIDATE_INPUT()
	mode = VALIDATE_INPUT.out.mode
	sample_fastq_ch = VALIDATE_INPUT.out.samples	
	
	// QC raw reads
	fastqc_ch = sample_fastq_ch.map { sample_id, fastqs -> 
										tuple(sample_id, fastqs, "raw") 
									} 
	FASTQC_RAW(fastqc_ch)		
	
	// Index reference
	//PREP_REFERENCE()
	star_index_ch = PREP_REFERENCE.out().star_index
	salmon_index_ch = PREP_REFERENCE.out().salmon_index
	
	// Quantify reads using SALMON
	SALMON_QUANT(sample_fastq_ch, salmon_index_ch)
	
	// Align reads using STAR
	STAR_ALIGN(sample_fastq_ch, star_index_ch)
	sample_bam_ch = STAR_ALIGN.out.bam_indexed
	
	// QC STAR alignments	
	RSEQC(sample_bam_ch, mode)	
	
	// Combine reports using MultiQC
	multiqc_ch = Channel.empty()
		.mix(FASTQC_RAW.out.fastqc_zip)             		// MultiQC parses the ZIPs
		.mix(SALMON_QUANT.out.salmon_quant.map { it[1] }) 	// The folder		
		.mix(STAR_ALIGN.out.gene_counts)        			// ReadsPerGene.out.tab
		.mix(STAR_ALIGN.out.sj_tab)             			// SJ.out.tab
		.mix(STAR_ALIGN.out.star_log)           			// Log.final.out		
		.mix(RSEQC.out.rseqc_logs)              			// The *.txt, *.log, *.r files
		.collect()	
	MULTIQC(multiqc_ch)
	
}






























/*
================================================================================
🌐 Environment Variable Handling: Bash/SGE vs Nextflow
================================================================================

1) Classic Bash + SGE (qsub) Behavior:

   - Exported variables in your shell are NOT automatically available inside qsub jobs.
     Example:
        export FILE_FORMAT="fq.gz"
        qsub my_job.sh   # FILE_FORMAT NOT automatically seen in job
   - To pass variables, you must use:
        qsub -v FILE_FORMAT="$FILE_FORMAT" my_job.sh
     or use `-V` to export the full environment.
   - Each job runs in a separate shell, so exported variables do not propagate automatically.

2) Nextflow Behavior:

   - Use `params` in `main.nf` or `nextflow.config` for global pipeline parameters.
     Example:
        params.file_format = "fq.gz"
   - Nextflow automatically injects `params` and channel values into all processes,
     even when using SGE, Slurm, or other executors.
   - No need for `export` or `-v`; processes can access the value directly:
        echo "File format is ${params.file_format}"
   - Channels and process inputs/outputs are also passed automatically.
   
3) Recommendation:

   - Define pipeline-wide variables in `params` at the top of `main.nf` or in config.
   - Update `params` dynamically during file detection or pipeline initialization.
   - This ensures all processes, whether single-end or paired-end, have access to
     consistent parameters, without worrying about environment propagation.

--- FASTQ file channel ---


--- Regex and file validation notes ---

 1. `~/.../` syntax creates a Groovy Pattern object (like `re.compile()` in Python)
    Example: VALID_PATTERN = ~/.*((_Tumor|_Normal))?.*((_R|_r)[12]).*\.f(q|astq)\.gz/

 2. `==~` operator matches the **entire string** against the regex
    Example: file.getName() ==~ VALID_PATTERN
    - Returns true only if the whole filename matches the pattern

 3. `=~` operator matches **any substring** against the regex
    Example: file.getName() =~ /_R1/
    - Returns a Matcher object if any part of the filename matches the pattern

 4. For full filename validation, use `==~ VALID_PATTERN`
    For classification by partial patterns (like read suffix _R1/_r1), use `=~ /pattern/`

 5. This is why we can safely separate R1/R2 and _r1/_r2 files without mixing uppercase/lowercase
    while validating the entire filename at the same time.



# Required by 01_fasterq_dump.sh
params.srr_file="${params.proj_dir}/SRR_Acc_List.txt"
params.expt="RNASeq"

# Required by 02_rename_fastqs.sh
params.map_file="${params.proj_dir}/rename_map.txt"