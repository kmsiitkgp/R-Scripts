#!/usr/bin/env nextflow
nextflow.enable.dsl=2 

include { VALIDATE_INPUT } 	from './modules/validate.nf'
include { FASTQC } 			from './modules/fastqc.nf'
include { PREP_REFERENCE } 	from './modules/prep_reference.nf'
include { SALMON_QUANT } 	from './modules/salmon_quant.nf'
include { STAR_ALIGN } 		from './modules/star_align.nf'
include { RSEQC } 			from './modules/rseqc.nf'
include { MULTIQC } 		from './modules/multiqc.nf'

workflow {
	
	//FASTERQ_DUMP()
	//RENAME FASTQS()
	
	// Validate fastq filenames
    VALIDATE_INPUT()
	mode = VALIDATE_INPUT.out.mode
	samples_ch = VALIDATE_INPUT.out.samples	
	
	// QC raw reads
	FASTQC(samples_ch)		
	
	// Index reference
	//PREP_REFERENCE()
	
	// Quantify reads using SALMON
	SALMON_QUANT(samples_ch)
	
	// Align reads using STAR
	STAR_ALIGN(samples_ch)
	bam_indexed_ch = STAR_ALIGN.out.bam_indexed
	
	// QC STAR alignments	
	RSEQC(bam_indexed_ch, mode)	
	
	// Combine reports using MultiQC
	multiqc_ch = Channel.empty()
		.mix(FASTQC.out.fastqc_zip)             			// MultiQC parses the ZIPs
		.mix(STAR_ALIGN.out.star_log)           			// Log.final.out
		.mix(STAR_ALIGN.out.sj_tab)             			// SJ.out.tab
		.mix(STAR_ALIGN.out.gene_counts)        			// ReadsPerGene.out.tab
		.mix(SALMON_QUANT.out.salmon_quant.map { it[1] }) 	// The folder
		.mix(RSEQC.out.rseqc_logs)              			// The *.txt, *.log, *.r files
		.collect()	
	MULTIQC(multiqc_ch)
	
}



/*
=========================================================================================
 NOTE: This script implicitly inherits all variables defined in the 'nextflow.config' 
 file located in this directory. Parameters prefixed with 'params.' do not need 
 to be declared locally here.

 PIPELINE LOGIC RULES:
 1. CHANNEL CREATION: Channels are created HERE, using the strings defined in config.
    (e.g., read_ch = Channel.fromPath(params.raw_fastq_pattern))
 2. IMPLICIT LOOPING: Nextflow will automatically 'loop' over the contents of a channel.
    You do not need 'for' loops; calling a Process on a Channel triggers the 'Batch' logic.
 3. PARALLELISM: Every item in a Channel is processed as an independent HPC job.

 CHANNEL CREATION LOGIC:
 1. Channel.fromPath(...) creates a Nextflow channel that emits each file individually.
    Suitable for streaming into processes for parallel execution.
 2. .collect() converts the stream into a single Channel emitting a single List 
    containing all files. This is used here for the validation block below.
 3. To feed processes later, we transform this data back into a structured channel.

 CODE BLOCKS : script vs shell
1) script: block with triple double quotes """  (DEFAULT)
 - Nextflow interpolates variables using ${}.
 - Bash variables MUST be escaped with \${} to prevent Nextflow from trying to evaluate them.
This is ideal for Nextflow scripts.

2) shell: block with triple single quotes '''  (COPY-PASTE FRIENDLY)
- Nextflow interpolates variables using !{} instead of ${}.
- Bash variables ${} work normally (no escaping needed).
This is ideal for copy-pasting existing Bash scripts.

script:
"""
echo "Sample ID is ${sample_id}"
CURRENT_DIR=\$PWD
echo "I am running in \$CURRENT_DIR"
"""

shell:
'''
echo "Sample ID is !{sample_id}"

# Bash variable (standard Bash syntax)
CURRENT_DIR=$PWD
echo "I am running in $CURRENT_DIR"
'''

INPUT types:
input:
    path index                 // A single/multiple reference file or directory
    val  mode                  // A simple scalar value (string, int, etc.)
    tuple val(id), path(fq)    // A sample "package": sample ID + FASTQ file
    env  TOOL_OPTS             // An environment variable injected into the task
	
1) Nextflow input example: val vs env.

Say you have sub_script.sh as below
#!/usr/bin/env bash
echo "Mode is: $mode"

This process will error because mode is defined as 'val' and not visibile for sub_script.sh	
process EXAMPLE_VAL {
    input:
        val mode

    shell:
    '''
    echo "Mode is !{mode}"
    bash sub_script.sh
    '''
}

This process will work because mode is defined as 'env' and visibile for sub_script.sh
process EXAMPLE_ENV {
    input:
        env MODE

    shell:
    '''
    echo "Mode is $MODE"
    bash sub_script.sh
    '''
}

2) Nextflow input example: val vs tuple

Suppose we have a sample metadata item:  ['single-end', 'R1', 'fastq.gz']
	
Channel.from( tuple('single-end','R1','fastq.gz') ) \
    .set { tuple_ch }    // Channel named tuple_ch emitting a tuple

process EXAMPLE_TUPLE {

    input:
        tuple val(layout), val(read), val(ext) from tuple_ch

    script:
    """
    echo "Using tuple inputs — separate values:"
    echo "Layout: ${layout}"
    echo "Read:   ${read}"
    echo "Suffix: ${ext}"
    """
}
Notes:
  - Each value has a clear name
  - Order is explicit
  - Safer and more maintainable
  - Recommended for structured metadata

Channel.from( ['single-end', 'R1', 'fastq.gz'] ) \
    .map { list -> list } \
    .set { val_ch }      // Channel named val_ch emitting one list

process EXAMPLE_VAL {

    input:
        val res from val_ch

    script:
    """
    echo "Using val(res) — single list:"
    echo "Layout: !{res[0]}"
    echo "Read:   !{res[1]}"
    echo "Suffix: !{res[2]}"
    """
}
Notes:
  - Works, but depends on list order
  - Harder to read and maintain
  - Good only for homogeneous or small groups
  
=========================================================================================
*/


















































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