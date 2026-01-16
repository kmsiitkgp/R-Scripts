// =========================================================================================
// GLOBAL PARAMETERS
// All variables defined in this 'params' block are automatically injected into 'main.nf'.
// They can be overridden via the command line (e.g., --proj_dir /new/path).

// CONFIGURATION DESIGN RULES:
// 1. DATA TYPE: Only define static strings, numbers, or booleans here (e.g., file paths).
// 2. NO CHANNELS: Do NOT try to create Channels (Channel.fromPath) in this file. 
//    Nextflow configuration is parsed before the execution engine starts.
// 3. OVERRIDES: Any value in this 'params' block can be overridden via terminal 
//    arguments (e.g., --raw_fastq_dir /new/path).
// =========================================================================================

// params.project is passed from command line `nextflow run main.nf --project Manish_22RV1_ARCaPM`

params {
	
	project 	= "test"
	expt 		= "RNASeq"
	proj_dir	= "${System.getenv('HOME')}/scratch/${expt}/${project}/"	
	assay		= "RNA"						// Required by Rscripts
	require_bam	= "false"					// Required by 02a_cellranger_count.sh
	
	// ==============================================================================
	// ⚙️ Define Project Specific Parameters 
	// ==============================================================================
	
	if (params.project == "Manish_22RV1_ARCaPM") {
		nsamples = 24
		species	 = "Homo sapiens"

	} else if (params.project == "Manish_22RV1_Xenograft") {   
		nsamples = 12
		species	 = "Homo sapiens"

	} else if (params.project == "Manish_22RV1_Xenograft2") {    
		nsamples = 8
		species	 = "Homo sapiens"

	} else if (params.project == "Vera") {   
		nsamples = 12
		species	 = "Homo sapiens"

	} else if (params.project == "Sandrine_Quad") {   
		nsamples = 15
		species	 = "Mus musculus"

	} else if (params.project == "Supriya_Quad") {   
		nsamples = 8
		species  = "Mus musculus"

	} else if (params.project == "Xinyi") {    
		nsamples = 48
		species	 = "Homo sapiens"

	} else {
		error "Unknown project: ${params.project}"
	}
	
	// Re-define species for bash without spaces
	if (species == "Homo sapiens") {
		species = "Human"
	} else if (species == "Mus musculus") {
		species = "Mouse"
	} else {
		error "❌ Undefined species: ${species}"
	}

	// ==============================================================================
	// ⚙️ Define Global Project Parameters
	// ==============================================================================
		
    //Find absolute path for conda env using `conda info --envs`
    my_env = '/home/user/anaconda3/envs/RNASEQ'	
    
	STAR_ARGS = [
		"--runMode", "alignReads",
		"--twopassMode", "Basic",
		"--quantMode", "GeneCounts",
		"--sjdbOverhang", "100",
		"--readFilesCommand", "zcat",
		"--outFilterMultimapNmax", "10",
		"--outFilterMismatchNmax", "999",
		"--outFilterMismatchNoverReadLmax", "0.04",
		"--outFilterType", "BySJout",
		"--alignEndsType", "Local",
		"--alignIntronMin", "20",
		"--alignIntronMax", "1000000",
		"--alignMatesGapMax", "1000000",
		"--outSAMunmapped", "Within",
		"--outSAMtype", "BAM", "SortedByCoordinate"
	]
	
	SALMON_ARGS = [
		"--libType", "A",
		"--gcBias",
		"--seqBias",
		"--posBias"
	]
	
	// ==============================================================================
	// 📁 Define Pipeline Directory Structure
	// ==============================================================================

	// --- Reference directories ---
	ref_dir   			= "${System.getenv('HOME')}/NGSTools/Reference_Genomes/${species}"				// directory with fa and gtf files of reference genome
	star_index_dir		= "${System.getenv('HOME')}/NGSTools/Reference_Genomes_STAR/${species}"			// directory with STAR indexed genome
	salmon_index_dir	= "${System.getenv('HOME')}/NGSTools/Reference_Genomes_SALMON/${species}"		// directory with SALMON indexed genome
	rseqc_index_dir		= "${System.getenv('HOME')}/NGSTools/Reference_Genomes_RSEQC/${species}"		// directory with RSEQC genome BED
	
	// --- Reference files ---
	// Grab the first FASTA file
	def refFastaFiles = file("${ref_dir}").listFiles().findAll { it.name.endsWith('.fa') }
	ref_fasta  = refFastaFiles[0].toString()

	// Grab the first GTF file
	def refGtfFiles = file("${ref_dir}").listFiles().findAll { it.name.endsWith('.gtf') }
	ref_gtf  = refGtfFiles[0].toString()

	gtf_basename = "${ref_gtf.tokenize('/')[-1].replace('.gtf','')}"   					// use tokenize('/')[-1] instead of basename
	salmon_transcritome_fasta 	= "${salmon_index_dir}/${gtf_basename}.transcripts.fa"	
	salmon_gentrome_fasta     	= "${salmon_index_dir}/${gtf_basename}.gentrome.fa"
	salmon_decoy              	= "${salmon_index_dir}/${gtf_basename}.decoy.txt"
	rseqc_bed					= "${rseqc_index_dir}/${gtf_basename}.bed"  

	// --- Result directories ---
	scripts_dir			= "${System.getenv('HOME')}/projects/${expt}/scripts/"			// directory with pipeline scripts

	log_dir				= "${proj_dir}/01.Logs/"

	fastq_dir			= "${proj_dir}/02.FastQ/"
	raw_fastq_dir		= "${fastq_dir}/1_raw/"				// directory with raw reads
	fastp_dir			= "${fastq_dir}/2_fastp/"        	// directory with fastp trimmed reads
	trimgalore_dir		= "${fastq_dir}/3_trimgalore/" 		// Keep for flexibility
	cutadapt_dir		= "${fastq_dir}/4_cutadapt/"    	// Keep for flexibility

	fastqc_dir			= "${proj_dir}/03.FastQC/"
	raw_qc_dir			= "${fastqc_dir}/1_raw/"			// directory with FastQC results of raw reads
	fastp_qc_dir		= "${fastqc_dir}/2_fastp/"			// directory with FastQC results of fastp trimmed reads
	trimgalore_qc_dir	= "${fastqc_dir}/3_trimgalore/"
	cutadapt_qc_dir		= "${fastqc_dir}/4_cutadapt/"

	star_dir			= "${proj_dir}/04.STAR/"
	rseqc_dir			= "${proj_dir}/05.RSeQC/"
	salmon_dir			= "${proj_dir}/06.Salmon/"
	deseq2_dir			= "${proj_dir}/07.DESeq2/"
	multiqc_dir			= "${proj_dir}/08.MultiQC/"
	
	multiqc_titlename	= "${params.project} MultiQC Report"
	multiqc_filename    = "${params.project}_MultiQC_Report"
	
}

// 2. Define the process configuration second (The "Instructions")
process {
    
	// Default settings for all processes
	// Only retry if it's a memory-related exit code (137, 139, 140, etc.)	
	errorStrategy 	= { task.exitStatus in [137, 139, 140, 143] ? 'retry' : 'terminate' }
    maxRetries    	= 3
	conda 		  	= "${params.my_env}"
	
	// Resource buckets
	withLabel: 'process_high' {
        cpus   = 8
        memory = { 48.GB * task.attempt }
    }

    withLabel: 'process_medium' {
        cpus   = 4
        memory = { 12.GB * task.attempt }
    }
	
	withLabel: 'process_low' {
        cpus   = 1
        memory = { 2.GB * task.attempt }
    }
	
	// SALMON settings
	withName: 'SALMON_QUANT' {		
		publishDir = [
            [ path: { "${params.salmon_dir}" },  			mode: 'copy' ], 
			[ path: { "${params.log_dir}" },     			mode: 'copy', pattern: "*.SALMON.error.log" ],		
			[ path: { "${params.salmon_dir}/quant_files" }, mode: 'copy', pattern: "**/*.sf",
			  saveAs: { filename -> 
			  // Add sample_id prefix
			  def sample = filename.tokenize('/')[0]   // assumes folder name = sample_id
			  return "${sample}.quant.sf" } ]
		]
	}	

    // STAR settings
    withName: 'STAR_ALIGN' {        
		publishDir = [
            [ path: { "${params.star_dir}" },  	mode: 'copy', 	pattern: "*.bam" ],
			[ path: { "${params.star_dir}" },  	mode: 'copy', 	pattern: "*.bai" ],
			[ path: { "${params.star_dir}" },  	mode: 'copy', 	pattern: "*.ReadsPerGene.out.tab" ],
			[ path: { "${params.star_dir}" },  	mode: 'copy', 	pattern: "*.Log.final.out" ],            
            [ path: { "${params.log_dir}" },   	mode: 'copy', 	pattern: "*.STAR.error.log" ]
		]     
    }
	
	// RSEQC settings
	withName: 'RSEQC' {    
		publishDir = [            
            [ path: params.rseqc_dir, 	mode: 'copy', 	pattern: "*.{pdf,jpeg,png,tiff}" ],
			[ path: params.rseqc_dir, 	mode: 'copy', 	pattern: "*.{txt,log,r}" ],
            [ path: params.log_dir,   	mode: 'copy', 	pattern: "*.RSEQC.error.log" ]
		]     
    }
	
	// MULTIQC settings
	withName: 'MULTIQC' {    
		publishDir = [            
            [ path: params.multiqc_dir,		mode: 'copy', 	pattern: "${params.multiqc_filename}.html" ],
			[ path: params.multiqc_dir,		mode: 'copy', 	pattern: "${params.multiqc_filename}_data" ],
            [ path: params.log_dir,			mode: 'copy', 	pattern: "MULTIQC.error.log" ]
		]     
    }
	
	
}

/*
executor {
    name = 'local'   // or 'slurm', 'pbs', 'sge', 'pbspro', 'awsbatch', etc.
    queueSize = 100  // max parallel jobs
}

conda {
    enabled = true
    autoMounts = true   // if using Singularity/remote paths
}
*/