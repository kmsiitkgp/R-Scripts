# 🚀 Complete Nextflow DSL2 Cheatsheet for Bioinformatics

> **Purpose**: Nextflow orchestrates workflows by creating channels and processes  
> **Key Concept**: Channels are asynchronous queues (pipes) that connect processes

---

## 📋 Table of Contents
1. [Nextflow Fundamentals](#1-nextflow-fundamentals)
2. [Channels (Data Pipes)](#2-channels-data-pipes)
3. [Processes (Compute Units)](#3-processes-compute-units)
4. [Workflows (Orchestration)](#4-workflows-orchestration)
5. [Process Inputs & Outputs](#5-process-inputs--outputs)
6. [PublishDir (Output Management)](#6-publishdir-output-management)
7. [Process Directives](#7-process-directives)
8. [Configuration (nextflow.config)](#8-configuration-nextflowconfig)
9. [Error Handling & Retry](#9-error-handling--retry)
10. [Executors (Local, SGE, AWS)](#10-executors-local-sge-aws)
11. [Containers (Docker/Singularity)](#11-containers-dockersingularity)
12. [Module System (Include/Import)](#12-module-system-includeimport)

---

## 1. Nextflow Fundamentals

### 1.1 DSL2 Syntax
```groovy
// Enable DSL2 (modern Nextflow syntax)
nextflow.enable.dsl=2

// DSL2 Features:
// - Explicit workflow blocks
// - Modular processes (one per file)
// - Reusable workflows
// - Named outputs (emit:)
```

### 1.2 Pipeline Structure
```bash
my_pipeline/
├── main.nf              # Main workflow entry point
├── nextflow.config      # Configuration file
├── modules/             # Process modules
│   ├── fastqc.nf
│   ├── star_align.nf
│   └── salmon_quant.nf
└── work/                # Execution directory (auto-created)
```

### 1.3 Running Pipelines
```bash
# Basic run
nextflow run main.nf

# With profile
nextflow run main.nf -profile sge,xinyi

# With resume (use cached results, workflow stops if terminal is closed)
nextflow run main.nf -profile xinyi -resume

# With resume (use cached results, workflow doesnt stop if terminal is closed)
nextflow run main.nf -profile xinyi -resume -bg

# Override parameters
nextflow run main.nf --project "NewProject" --cpus 16

# Generate reports
nextflow run main.nf -with-report report.html -with-timeline timeline.html -with-dag dag.png

# Clean work directory
nextflow clean -f
```

---

## 2. Channels (Data Pipes)

### 2.1 Channel Factories

#### 2.1.1 Channel.fromPath (File Streaming)
```groovy
// Create channel from file glob pattern (LAZY - creates stream)
fastq_ch = Channel.fromPath("*.fastq.gz")
// Emits: file1.fastq.gz, file2.fastq.gz, file3.fastq.gz (one at a time)

// With options
fastq_ch = Channel.fromPath("*.fastq.gz", checkIfExists: true)

// Recursive search
all_fastq = Channel.fromPath("**/*.fastq.gz")  // Searches subdirectories
```

**Real Example from validate.nf:**
```groovy
// GROOVY: Immediate file collection (returns Set)
fastq_files_set = files("${fastq_dir}/*.f*q.gz")  
fastq_files_list = fastq_files_set.collect()      // Convert to List
fastq_files = fastq_files_list.sort()             // Sort

// NEXTFLOW: Lazy channel creation
fastq_ch = Channel.fromPath("${fastq_dir}/*.f*q.gz")
```

#### 2.1.2 Channel.fromFilePairs (Auto-pairing R1/R2)
```groovy
// Automatically pair R1/R2 files
paired_ch = Channel.fromFilePairs("*_{R1,R2}.fastq.gz")
// Emits: [sample1, [sample1_R1.fastq.gz, sample1_R2.fastq.gz]]
//        [sample2, [sample2_R1.fastq.gz, sample2_R2.fastq.gz]]

// With custom pattern
paired_ch = Channel.fromFilePairs("*_{1,2}.fq.gz", flat: true)
// flat: true -> [sample1, R1.fq.gz, R2.fq.gz] (no nested list)
// flat: false -> [sample1, [R1.fq.gz, R2.fq.gz]] (default)
```

#### 2.1.3 Channel.value (Singleton/Constant)
```groovy
// Value channel (emits infinitely, reusable)
ref_genome = Channel.value(file("genome.fa", checkIfExists: true))

// Used when one input is shared by many samples
// Example: All samples use the same reference genome
```

**Real Example from main.nf:**
```groovy
// Create value channels for reference files
ref_fasta_ch = Channel.value(file(params.ref_fasta[params.species], checkIfExists: true))
ref_gtf_ch   = Channel.value(file(params.ref_gtf[params.species], checkIfExists: true))

// These channels can be used by EVERY sample without being consumed
PREP_REFERENCE(ref_fasta_ch, ref_gtf_ch)
```

#### 2.1.4 Channel.of (Literal Values)
```groovy
// Create channel from literal values
numbers_ch = Channel.of(1, 2, 3, 4, 5)
// Emits: 1, 2, 3, 4, 5

samples_ch = Channel.of('Sample1', 'Sample2', 'Sample3')
// Emits: 'Sample1', 'Sample2', 'Sample3'
```

#### 2.1.5 Channel.fromList (From Groovy List)
```groovy
// Convert Groovy list to Nextflow channel
def my_list = ["A", "B", "C"]
list_ch = Channel.fromList(my_list)
// Emits: "A", "B", "C"
```

#### 2.1.6 Channel.empty (Placeholder)
```groovy
// Create empty channel (useful for conditional logic)
empty_ch = Channel.empty()

// Example: Mix outputs conditionally
multiqc_ch = Channel.empty()
    .mix(FASTQC.out.zip)
    .mix(STAR.out.log)
    .collect()
```

**Real Example from main.nf:**
```groovy
// Start with empty channel, then mix in all QC outputs
multiqc_ch = Channel.empty()
    .mix(FASTQC_RAW.out.fastqc_zip)                     // FastQC ZIPs
    .mix(SALMON_QUANT.out.salmon_quant.map { it[1] })   // SALMON dirs
    .mix(STAR_ALIGN.out.gene_counts)                    // Gene counts
    .mix(STAR_ALIGN.out.sj_tab)                         // Junctions
    .mix(STAR_ALIGN.out.star_log)                       // Logs
    .mix(RSEQC.out.rseqc_logs)                          // RSeQC data
    .collect()  // Wait for all, then emit single list

MULTIQC(multiqc_ch)
```

### 2.2 Channel Operators

#### 2.2.1 .map (Transform Items)
```groovy
// Transform each item in channel
fastq_ch = Channel.fromPath("*.fastq.gz")
    .map { file -> [file.simpleName, file] }
// Input:  sample1.fastq.gz
// Output: [sample1, sample1.fastq.gz]
```

**Real Example from main.nf:**
```groovy
// Add "raw" tag for FastQC output organization
fastqc_ch = sample_fastq_ch.map { sample_id, fastqs -> 
    tuple(sample_id, fastqs, "raw") 
}
// Input:  [Sample1, [R1.fq.gz, R2.fq.gz]]
// Output: [Sample1, [R1.fq.gz, R2.fq.gz], "raw"]

FASTQC_RAW(fastqc_ch)
```

**Real Example from validate.nf:**
```groovy
// Extract sample_id and pair R1/R2 files
grouped_samples_ch = R1_FASTQS_ch.map { r1 ->
    def idx = r1.name.lastIndexOf(Read1_TAG)
    def sample_id = (idx != -1) ? r1.name.take(idx) : r1.simpleName
    
    if (MODE == "PAIRED_END") {
        def r2_name = r1.name.reverse()
                             .replaceFirst(Read1_TAG.reverse(), Read2_TAG.reverse())
                             .reverse()
        def r2 = r1.parent.resolve(r2_name)
        return [sample_id, [r1, r2]]
    } else {
        return [sample_id, [r1]]
    }
}
```

#### 2.2.2 .collect (Wait for All, Emit List)
```groovy
// Collect all items into a single list
fastq_ch = Channel.fromPath("*.fastq")
    .collect()
// Input:  file1.fastq, file2.fastq, file3.fastq (streamed)
// Output: [file1.fastq, file2.fastq, file3.fastq] (one emission)
```

**Real Example from main.nf:**
```groovy
// Collect outputs to ensure processes wait for all samples
star_index_ch   = PREP_REFERENCE.out.star_index_dir.collect()
salmon_index_ch = PREP_REFERENCE.out.salmon_index_dir.collect()
ref_bed_ch      = PREP_REFERENCE.out.ref_bed.collect()

// Without .collect(), downstream processes might start before index is ready
```

#### 2.2.3 .flatten (Expand Lists)
```groovy
// Break nested lists into individual items
nested_ch = Channel.of([1, 2], [3, 4], [5, 6])
    .flatten()
// Input:  [1, 2], [3, 4], [5, 6]
// Output: 1, 2, 3, 4, 5, 6
```

#### 2.2.4 .filter (Conditional Selection)
```groovy
// Keep only items matching condition
numbers_ch = Channel.of(1, 2, 3, 4, 5, 6, 7, 8, 9, 10)
    .filter { it % 2 == 0 }
// Output: 2, 4, 6, 8, 10
```

#### 2.2.5 .unique (Remove Duplicates)
```groovy
// Remove duplicate items
samples_ch = Channel.of("A", "B", "A", "C", "B", "D")
    .unique()
// Output: "A", "B", "C", "D"
```

#### 2.2.6 .view (Debug Print)
```groovy
// Print items as they flow through channel (debugging)
fastq_ch = Channel.fromPath("*.fastq")
    .view()  // Prints each filename
    .map { file -> [file.simpleName, file] }
    .view { "Processed: $it" }  // Custom message
```

#### 2.2.7 .join (Merge by Key)
```groovy
// Join two channels by matching keys
ch1 = Channel.of(
    ['sample1', 'A'],
    ['sample2', 'B']
)

ch2 = Channel.of(
    ['sample1', 'X'],
    ['sample2', 'Y']
)

joined = ch1.join(ch2)
// Output: [sample1, A, X]
//         [sample2, B, Y]
```

#### 2.2.8 .combine (Cartesian Product)
```groovy
// Every item from ch1 with every item from ch2
samples = Channel.of('A', 'B')
tools = Channel.of('tool1', 'tool2')

combinations = samples.combine(tools)
// Output: [A, tool1]
//         [A, tool2]
//         [B, tool1]
//         [B, tool2]
```

#### 2.2.9 .groupTuple (Group by Key)
```groovy
// Group multiple values by the same key
ch = Channel.of(
    ['sample1', 'file1.bam'],
    ['sample1', 'file2.bam'],
    ['sample2', 'file3.bam']
)

grouped = ch.groupTuple()
// Output: [sample1, [file1.bam, file2.bam]]
//         [sample2, [file3.bam]]
```

#### 2.2.10 .mix (Merge Channels)
```groovy
// Combine multiple channels into one
ch1 = Channel.of(1, 2, 3)
ch2 = Channel.of(4, 5, 6)

mixed = ch1.mix(ch2)
// Output: 1, 2, 3, 4, 5, 6 (order not guaranteed)
```

---

## 3. Processes (Compute Units)

### 3.1 Basic Process Structure
```groovy
process PROCESS_NAME {
    // Directives (optional)
    tag "..."
    label 'process_low'
    
    // Inputs
    input:
    <input declarations>
    
    // Outputs
    output:
    <output declarations>
    
    // Execution script
    script:
    """
    bash commands here
    """
}
```

### 3.2 Complete Process Example

**Real Example from fastqc.nf:**
```groovy
process FASTQC {
    tag "FastQC on ${sample_id} (${read_type})"
    label 'process_low'
    
    input:
    tuple val(sample_id), path(fastq_files), val(read_type)
    
    output:
    path("*_fastqc.zip"),        emit: fastqc_zip
    path("*_fastqc.html"),       emit: fastqc_html
    path("*.FASTQC.error.log"),  emit: fastqc_error_log
    
    script:
    def LOG = "${sample_id}.${read_type}.FASTQC.error.log"
    
    """
    fastqc \
        --threads "${task.cpus}" \
        --quiet \
        ${fastq_files.join(' ')} \
        1>> "${LOG}" 2>&1 \
        || { echo "❌ ERROR: FastQC failed" | tee -a "${LOG}" >&2; exit 1; }
    
    echo "✅ SUCCESS: FastQC completed" >> "${LOG}"
    """
}
```

### 3.3 Script vs Shell Block

#### 3.3.1 script: Block (Default)
```groovy
process EXAMPLE_SCRIPT {
    input:
    val sample_id
    
    script:
    """
    # Nextflow variables use ${variable}
    echo "Sample ID is ${sample_id}"
    
    # Bash variables MUST be escaped with backslash
    CURRENT_DIR=\$PWD
    echo "Running in \$CURRENT_DIR"
    """
}
```

#### 3.3.2 shell: Block (Copy-Paste Friendly)
```groovy
process EXAMPLE_SHELL {
    input:
    val sample_id
    
    shell:
    '''
    # Nextflow variables use !{variable}
    echo "Sample ID is !{sample_id}"
    
    # Bash variables work normally (no escaping)
    CURRENT_DIR=$PWD
    echo "Running in $CURRENT_DIR"
    '''
}
```

**When to use which:**
- `script:` → Default, use when mixing Nextflow and Bash variables
- `shell:` → Use when copy-pasting existing Bash scripts

---

## 4. Workflows (Orchestration)

### 4.1 Basic Workflow Structure
```groovy
workflow {
    // Create channels
    fastq_ch = Channel.fromPath("*.fastq")
    
    // Call processes
    PROCESS_A(fastq_ch)
    
    // Use outputs
    PROCESS_B(PROCESS_A.out.result)
}
```

### 4.2 Named Workflows (Reusable)</p>
```groovy
workflow MY_WORKFLOW {
    take:
    input_ch    // Required input
    
    main:
    PROCESS_A(input_ch)
    PROCESS_B(PROCESS_A.out)
    
    emit:
    result = PROCESS_B.out.final
}

// Call from main workflow
workflow {
    data_ch = Channel.fromPath("*.txt")
    MY_WORKFLOW(data_ch)
}
```

**Real Example from validate.nf:**
```groovy
workflow VALIDATE_INPUT {
    take:
    data_dir  // Input: directory path
    
    main:
    // Validation logic here
    fastq_files = files("${data_dir}/*.f*q.gz")
    // ... validation code ...
    
    emit:
    mode             = MODE
    samples          = grouped_samples_ch
    total_samples    = TOTAL_SAMPLES
}

// Usage in main.nf
workflow {
    VALIDATE_INPUT(params.raw_fastq_dir())
    
    mode            = VALIDATE_INPUT.out.mode.collect()
    sample_fastq_ch = VALIDATE_INPUT.out.samples
}
```

---

## 5. Process Inputs & Outputs

### 5.1 Input Types

#### 5.1.1 val (Simple Values)
```groovy
process EXAMPLE {
    input:
    val sample_id      // String, integer, boolean, etc.
    val cpus           // Number
    
    script:
    """
    echo "Processing ${sample_id} with ${cpus} cores"
    """
}
```

#### 5.1.2 path (Files/Directories)
```groovy
process EXAMPLE {
    input:
    path fastq         // Single file
    path(fastqs)       // List of files
    path index_dir     // Directory
    
    script:
    """
    # Nextflow stages files in work directory
    ls -lh ${fastq}
    ls -lh ${index_dir}
    """
}
```

#### 5.1.3 tuple (Grouped Values)
```groovy
process EXAMPLE {
    input:
    tuple val(sample_id), path(fastq_files)
    // Receives: [sample1, [R1.fq.gz, R2.fq.gz]]
    
    script:
    """
    echo "Sample: ${sample_id}"
    echo "Files: ${fastq_files.join(', ')}"
    """
}
```

**Real Example from star_align.nf:**
```groovy
process STAR_ALIGN {
    input:
    tuple val(sample_id), path(fastq_files)  // [ID, [files]]
    path(star_index_dir)                     // Shared index
    
    script:
    def MATES_ARGS = fastq_files.size() == 2 ?
        "--readFilesIn ${fastq_files[0]} ${fastq_files[1]}" :
        "--readFilesIn ${fastq_files[0]}"
    
    """
    STAR --genomeDir ${star_index_dir} ${MATES_ARGS} ...
    """
}
```

#### 5.1.4 env (Environment Variables)
```groovy
process EXAMPLE {
    input:
    env MODE           // Set as environment variable
    
    script:
    """
    # MODE is available to all subprocesses
    echo "Mode is: $MODE"
    bash sub_script.sh  # Can access $MODE
    """
}
```

**val vs env:**
- `val` → Only visible in Nextflow interpolation `${var}`
- `env` → Exported to shell, visible in subscripts

### 5.2 Output Types

#### 5.2.1 path (Files/Directories)
```groovy
process EXAMPLE {
    output:
    path "*.bam"                    // Glob pattern
    path "${sample_id}.txt"         // Specific file
    path "output_dir"               // Directory
    
    script:
    """
    touch file1.bam file2.bam
    echo "data" > ${sample_id}.txt
    mkdir output_dir
    """
}
```

#### 5.2.2 Named Outputs (emit:)
```groovy
process STAR_ALIGN {
    output:
    path "*.bam",         emit: bam_file
    path "*.bai",         emit: bam_index
    path "*.tab",         emit: gene_counts
    
    script:
    """
    # Generate outputs
    """
}

// Usage in workflow
workflow {
    STAR_ALIGN(inputs)
    
    // Access named outputs
    bam_ch = STAR_ALIGN.out.bam_file
    idx_ch = STAR_ALIGN.out.bam_index
}
```

**Real Example from star_align.nf:**
```groovy
process STAR_ALIGN {
    output:
    tuple val(sample_id),
        path("${sample_id}.Aligned.sortedByCoord.out.bam"),
        path("${sample_id}.Aligned.sortedByCoord.out.bam.bai"),
        emit: bam_indexed  // Tuple output
    
    path "${sample_id}.ReadsPerGene.out.tab",   emit: gene_counts
    path "${sample_id}.SJ.out.tab",             emit: sj_tab
    path "${sample_id}.Log.final.out",          emit: star_log
}

// Usage
workflow {
    STAR_ALIGN(fastq_ch, index_ch)
    
    sample_bam_ch = STAR_ALIGN.out.bam_indexed  // [id, bam, bai]
    counts_ch     = STAR_ALIGN.out.gene_counts
}
```

#### 5.2.3 optional (Conditional Outputs)
```groovy
process EXAMPLE {
    output:
    path "*.pdf",  optional: true   // May not exist
    
    script:
    """
    # Only create PDF for certain conditions
    if [ condition ]; then
        generate_plot.py > plot.pdf
    fi
    """
}
```

**Real Example from rseqc.nf:**
```groovy
process RSEQC {
    output:
    path("${sample_id}*.{pdf,jpeg,png}"),  emit: rseqc_plots,  optional: true
    path("${sample_id}*.{txt,log,r}"),     emit: rseqc_logs,   optional: true
    
    script:
    """
    # Some plots only generated for paired-end data
    if [[ "${mode}" == "PAIRED_END" ]]; then
        inner_distance.py ... > plot.pdf
    fi
    """
}
```

---

## 6. PublishDir (Output Management)

### 6.1 Basic PublishDir
```groovy
process EXAMPLE {
    publishDir "/results", mode: 'copy'
    
    output:
    path "*.txt"
    
    script:
    """
    echo "data" > result.txt
    """
}
```

### 6.2 PublishDir Modes

```groovy
// COPY - Safe, uses disk space
publishDir "/results", mode: 'copy'

// SYMLINK - Saves space, fragile
publishDir "/results", mode: 'symlink'

// MOVE - Fast, removes from work/
publishDir "/results", mode: 'move'

// LINK - Hard link (same filesystem only)
publishDir "/results", mode: 'link'
```

### 6.3 publishDir in nextflow.config (RECOMMENDED)

**Real Example from nextflow.config:**
```groovy
process {
    withName: 'FASTQC' {
        publishDir = [
            // HTML reports
            [
                path: { "${params.proj_dir()}/02.FastQC/${read_type}" },
                mode: 'copy',
                pattern: "*.html"
            ],
            // ZIP files for MultiQC
            [
                path: { "${params.proj_dir()}/02.FastQC/${read_type}" },
                mode: 'copy',
                pattern: "*.zip"
            ],
            // Error logs
            [
                path: "${params.proj_dir()}/07.Logs/",
                mode: 'copy',
                pattern: "*.FASTQC.error.log"
            ]
        ]
    }
}
```

### 6.4 Dynamic publishDir Paths

**Real Example - SALMON_QUANT:**
```groovy
withName: 'SALMON_QUANT' {
    publishDir = [
        // Full output directory
        [
            path: "${params.proj_dir()}/03.Salmon/",
            mode: 'copy',
            pattern: { "${sample_id}" }  // Closure with sample_id
        ],
        // Just quant.sf, renamed
        [
            path: "${params.proj_dir()}/03.Salmon/quant_files",
            mode: 'copy',
            pattern: "**/quant.sf",
            saveAs: { filename -> "${sample_id}.quant.sf" }
        ]
    ]
}
```

### 6.5 saveAs (Rename on Publish)
```groovy
publishDir "/results", saveAs: { filename ->
    filename.endsWith('.bam') ? "bam_files/${filename}" : filename
}

// Or in config
saveAs: { filename -> "${sample_id}.${filename}" }
```

---

## 7. Process Directives

### 7.1 Resource Allocation

```groovy
process EXAMPLE {
    cpus 8              // Number of CPU cores
    memory '32 GB'      // Memory allocation
    time '2h'           // Maximum runtime
    disk '100 GB'       // Disk space
    
    script:
    """
    # Use task variables
    my_tool --threads ${task.cpus} --memory ${task.memory.toGiga()}G
    """
}
```

### 7.2 Labels (Group Resources)

**In process:**
```groovy
process FASTQC {
    label 'process_low'    // Reference to config
    
    script:
    """
    fastqc --threads ${task.cpus}
    """
}
```

**In nextflow.config:**
```groovy
process {
    withLabel: 'process_low' {
        cpus   = 2
        memory = { 3.GB * task.attempt }
    }
    
    withLabel: 'process_medium' {
        cpus   = 4
        memory = { 12.GB * task.attempt }
    }
    
    withLabel: 'process_high' {
        cpus   = 8
        memory = { 48.GB * task.attempt }
    }
}
```

### 7.3 tag (Display in Logs)
```groovy
process STAR_ALIGN {
    tag "$sample_id"    // Shows in Nextflow logs
    
    input:
    tuple val(sample_id), path(fastq)
    
    script:
    """
    STAR ...
    """
}

// Log output:
// [STAR_ALIGN] Sample1
// [STAR_ALIGN] Sample2
```

### 7.4 container (Singularity/Docker)
```groovy
process FASTQC {
    container 'quay.io/biocontainers/fastqc:0.12.1--hdfd78af_0'
    
    script:
    """
    fastqc ...
    """
}
```

### 7.5 cache (Control Caching)
```groovy
process EXAMPLE {
    cache 'deep'      // Cache based on file content (default: 'standard')
    cache false       // Disable caching
    
    script:
    """
    ...
    """
}
```

### 7.6 errorStrategy (Error Handling)
```groovy
process EXAMPLE {
    errorStrategy 'retry'      // Retry on failure
    errorStrategy 'ignore'     // Continue despite failure
    errorStrategy 'terminate'  // Stop pipeline (default)
    maxRetries 3               // Max retry attempts
    
    script:
    """
    ...
    """
}
```

**Real Example from nextflow.config:**
```groovy
process {
    // Retry only on memory errors
    errorStrategy = { task.exitStatus in [137, 139, 140, 143] ? 'retry' : 'terminate' }
    maxRetries = 3
}
```

### 7.7 when (Conditional Execution)
```groovy
process OPTIONAL_STEP {
    when:
    params.run_optional == true
    
    input:
    path data
    
    script:
    """
    # Only runs if when: condition is true
    """
}
```

---

## 8. Configuration (nextflow.config)

### 8.1 params (Pipeline Parameters)
```groovy
params {
    project = "MyProject"
    species = "Human"
    cpus    = 8
    
    // Closures for dynamic evaluation
    proj_dir = { "${params.base_dir}/${params.project}" }
}

// Access in workflow
println params.project       // "MyProject"
println params.proj_dir()    // Evaluates closure
```

### 8.2 profiles (Environment Configs)
```groovy
profiles {
    standard {
        process.executor = 'local'
    }
    
    sge {
        process.executor = 'sge'
        process.queue    = 'all.q'
        process.penv     = 'smp'
        workDir = "/scratch/work"
        
        executor {
            queueSize = 10    // Max concurrent jobs
        }
    }
    
    aws {
        process.executor = 'awsbatch'
        process.queue    = 'my-queue'
        workDir = 's3://my-bucket/work'
    }
}

// Usage: nextflow run main.nf -profile sge
```

**Real Example:**
```groovy
profiles {
    xinyi {
        params.project        = "Xinyi"
        params.species        = "Human"
        params.genome_version = "GRCh38"
    }
    
    sandrine {
        params.project        = "Sandrine_Quad"
        params.species        = "Mouse"
        params.genome_version = "GRCm39"
    }
}
```

### 8.3 process Configuration
```groovy
process {
    // Default for all processes
    cpus   = 1
    memory = 2.GB
    
    // Process-specific
    withName: 'STAR_ALIGN' {
        cpus   = 8
        memory = '48 GB'
        container = 'quay.io/biocontainers/star:2.7.11b'
    }
    
    // Label-based
    withLabel: 'process_high' {
        cpus   = 8
        memory = { 48.GB * task.attempt }
    }
}
```

### 8.4 Singularity Configuration
```groovy
singularity {
    enabled     = true
    autoMounts  = false
    runOptions  = "--bind /scratch --bind /home/user"
    cacheDir    = "/scratch/singularity_cache"
}
```

### 8.5 Trace, Report, Timeline
```groovy
trace {
    enabled   = true
    file      = 'trace.txt'
    overwrite = true
}

report {
    enabled   = true
    file      = 'report.html'
    overwrite = true
}

timeline {
    enabled   = true
    file      = 'timeline.html'
    overwrite = true
}
```

---

## 9. Error Handling & Retry

### 9.1 Basic Error Handling
```groovy
process EXAMPLE {
    errorStrategy 'retry'
    maxRetries 3
    
    script:
    """
    command_that_might_fail \
        || { echo "Error occurred"; exit 1; }
    """
}
```

### 9.2 Memory-Based Retry Strategy

**Real Example from nextflow.config:**
```groovy
process {
    // Only retry on memory errors
    errorStrategy = { task.exitStatus in [137, 139, 140, 143] ? 'retry' : 'terminate' }
    maxRetries = 3
    
    // Memory scales with attempt
    withLabel: 'process_high' {
        memory = { 48.GB * task.attempt }
        // Attempt 1: 48GB
        // Attempt 2: 96GB
        // Attempt 3: 144GB
    }
}
```

### 9.3 Exit Code Reference
```groovy
// Common exit codes
// 137 = SIGKILL (out of memory)
// 139 = SIGSEGV (segmentation fault)
// 140 = SIGTERM with stack overflow
// 143 = SIGTERM (terminated, often memory)
```

### 9.4 Error Logging Pattern

**Real Example from processes:**
```groovy
script:
def LOG = "${sample_id}.PROCESS.error.log"

"""
tool --input ${input} \
    1>> "${LOG}" 2>&1 \
    || { echo "❌ ERROR: Tool failed for ${sample_id}" | tee -a "${LOG}" >&2; exit 1; }

echo "✅ SUCCESS: Tool completed for ${sample_id}" >> "${LOG}"
"""
```

**Explanation:**
- `1>> "${LOG}"` → Redirect stdout to log (append)
- `2>&1` → Redirect stderr to stdout (both go to log)
- `|| { ... }` → Execute on non-zero exit code
- `tee -a "${LOG}" >&2` → Write to log AND stderr (visible in Nextflow)
- `exit 1` → Signal failure to Nextflow

---

## 10. Executors (Local, SGE, AWS)

### 10.1 Local Executor (Default)
```groovy
process.executor = 'local'

// Runs processes on the machine running Nextflow
// Good for: testing, small datasets, local workstations
```

### 10.2 SGE (Sun Grid Engine)

**Real Example from nextflow.config:**
```groovy
profiles {
    sge {
        process.executor       = 'sge'
        process.queue          = 'all.q'
        process.penv           = 'smp'
        process.clusterOptions = { "-S /bin/bash -l h_vmem=${task.memory.toGiga()}G" }
        workDir                = "${params.base_dir}/work"
        
        executor {
            queueSize = 10  // Max jobs in queue
        }
    }
}
```

**Usage:**
```bash
nextflow run main.nf -profile sge,xinyi -resume
```

### 10.3 SLURM Executor
```groovy
profiles {
    slurm {
        process.executor = 'slurm'
        process.queue    = 'normal'
        process.time     = '2h'
        process.clusterOptions = '--account=myproject'
        
        executor {
            queueSize = 50
        }
    }
}
```

### 10.4 AWS Batch

**Real Example from nextflow.config:**
```groovy
profiles {
    aws {
        process.executor = 'awsbatch'
        process.queue    = 'your-job-queue-name'
        process.memory   = { task.memory.toMega() }
        workDir          = "s3://my-nextflow-bucket/work"
        
        aws {
            region = 'us-east-1'
            batch {
                jobRole = 'arn:aws:iam::ACCOUNT_ID:role/BatchServiceRole'
                cliPath = '/home/ec2-user/miniconda/bin/aws'
                maxParallelTransfers = 5
            }
        }
        
        // AWS Batch requires Docker
        docker.enabled      = true
        singularity.enabled = false
    }
}
```

### 10.5 Executor Comparison

| Executor | Use Case | Pros | Cons |
|----------|----------|------|------|
| `local` | Testing, small jobs | Simple, no setup | Limited resources |
| `sge` | HPC clusters | Good for legacy systems | Complex setup |
| `slurm` | Modern HPC | Widely used, efficient | Requires SLURM cluster |
| `awsbatch` | Cloud computing | Scalable, elastic | Costs money, complex setup |
| `k8s` | Kubernetes | Modern, portable | Complex configuration |

---

## 11. Containers (Docker/Singularity)

### 11.1 Singularity (HPC Standard)

**In nextflow.config:**
```groovy
singularity {
    enabled     = true
    autoMounts  = false
    runOptions  = "--bind /scratch --bind /home/user"
    cacheDir    = "/scratch/singularity_cache"
}
```

**In process:**
```groovy
process FASTQC {
    container 'quay.io/biocontainers/fastqc:0.12.1--hdfd78af_0'
    
    script:
    """
    fastqc --version
    fastqc ${fastq}
    """
}
```

### 11.2 Docker (Local/Cloud)
```groovy
docker {
    enabled = true
    runOptions = '-u $(id -u):$(id -g)'
}

process EXAMPLE {
    container 'ubuntu:20.04'
    
    script:
    """
    echo "Running in Docker"
    """
}
```

### 11.3 Container Sources

```groovy
// BioContainers (recommended for bioinformatics)
container 'quay.io/biocontainers/fastqc:0.12.1--hdfd78af_0'
container 'quay.io/biocontainers/star:2.7.11b--h5ca1c30_4'
container 'quay.io/biocontainers/salmon:1.10.3--h7e5ed60_0'

// Docker Hub
container 'ubuntu:20.04'
container 'continuumio/miniconda3:latest'

// Local Singularity image
container '/path/to/image.sif'
```

### 11.4 Per-Process Containers

**Real Example from nextflow.config:**
```groovy
process {
    withName: 'FASTQC' {
        container = 'quay.io/biocontainers/fastqc:0.12.1--hdfd78af_0'
    }
    
    withName: 'STAR_ALIGN' {
        container = 'quay.io/biocontainers/star:2.7.11b--h5ca1c30_4'
    }
    
    withName: 'SALMON_QUANT' {
        container = 'quay.io/biocontainers/salmon:1.10.3--h7e5ed60_0'
    }
}
```

### 11.5 Bind Paths (Singularity)

**Common Issues:**
```groovy
// Problem: Container can't see files
singularity {
    autoMounts = true  // Let Nextflow guess (may fail)
}

// Solution: Explicit bind paths
singularity {
    autoMounts  = false
    runOptions  = "--bind /scratch --bind /home/user --bind /data"
}
```

**Troubleshooting Bind Paths:**
```bash
# Check symlink in work directory
ls -lh ~/scratch/work/3a/517fad.../.../file.fastq.gz
# lrwxrwxrwx ... -> /scratch/user/project/file.fastq.gz

# If Singularity sees different path:
# Container: /scratch/user/file.fastq.gz
# Actual:    /hpc/home/user/scratch/file.fastq.gz

# Solution: Bind both paths
runOptions = "--bind /scratch --bind /hpc/home/user"
```

---

## 12. Module System (Include/Import)

### 12.1 Basic Include
```groovy
// main.nf
include { FASTQC } from './modules/fastqc.nf'
include { STAR_ALIGN } from './modules/star_align.nf'

workflow {
    fastq_ch = Channel.fromPath("*.fastq")
    FASTQC(fastq_ch)
    STAR_ALIGN(FASTQC.out.trimmed)
}
```

### 12.2 Include with Alias
```groovy
// Use same process twice with different names
include { FASTQC as FASTQC_RAW } from './modules/fastqc.nf'
include { FASTQC as FASTQC_TRIMMED } from './modules/fastqc.nf'

workflow {
    raw_ch = Channel.fromPath("*_raw.fastq")
    trimmed_ch = Channel.fromPath("*_trimmed.fastq")
    
    FASTQC_RAW(raw_ch)
    FASTQC_TRIMMED(trimmed_ch)
}
```

**Real Example from main.nf:**
```groovy
include { FASTQC as FASTQC_RAW }     from './modules/fastqc.nf'
include { FASTQC as FASTQC_TRIMMED } from './modules/fastqc.nf'

workflow {
    // Raw reads QC
    fastqc_ch = sample_fastq_ch.map { id, fastqs -> tuple(id, fastqs, "raw") }
    FASTQC_RAW(fastqc_ch)
    
    // Trimmed reads QC (if trimming is done)
    // FASTQC_TRIMMED would be called here
}
```

### 12.3 Include Workflows
```groovy
// modules/validate.nf
workflow VALIDATE_INPUT {
    take:
    data_dir
    
    main:
    // Validation code
    
    emit:
    samples = sample_ch
    mode = detected_mode
}

// main.nf
include { VALIDATE_INPUT } from './modules/validate.nf'

workflow {
    VALIDATE_INPUT(params.raw_fastq_dir())
    
    samples_ch = VALIDATE_INPUT.out.samples
    mode_ch    = VALIDATE_INPUT.out.mode
}
```

### 12.4 Multiple Includes
```groovy
include { 
    FASTQC;
    MULTIQC;
    TRIMMOMATIC 
} from './modules/qc.nf'

include {
    STAR_ALIGN;
    SALMON_QUANT
} from './modules/align.nf'
```

---

## 13. Advanced Patterns

### 13.1 Conditional Process Execution

```groovy
workflow {
    fastq_ch = Channel.fromPath("*.fastq")
    
    // Option 1: when directive
    OPTIONAL_PROCESS(fastq_ch)
    
    // Option 2: Groovy if
    if (params.run_trimming) {
        TRIMMING(fastq_ch)
        next_ch = TRIMMING.out.trimmed
    } else {
        next_ch = fastq_ch
    }
    
    ALIGN(next_ch)
}

process OPTIONAL_PROCESS {
    when:
    params.run_optional
    
    script:
    """
    echo "Running optional step"
    """
}
```

### 13.2 Dynamic Input Creation

**Real Example from star_align.nf:**
```groovy
script:
// Create different arguments based on file count
def MATES_ARGS = fastq_files.size() == 2 ?
    "--readFilesIn ${fastq_files[0]} ${fastq_files[1]}" :
    "--readFilesIn ${fastq_files[0]}"

"""
STAR --genomeDir ${index} ${MATES_ARGS} --threads ${task.cpus}
"""
```

### 13.3 Multi-Path Publishing

**Real Example - Publishing to Multiple Locations:**
```groovy
withName: 'STAR_ALIGN' {
    publishDir = [
        // BAM files
        [
            path: "${params.proj_dir()}/04.STAR/",
            mode: 'copy',
            pattern: "*.bam"
        ],
        // BAM indexes
        [
            path: "${params.proj_dir()}/04.STAR/",
            mode: 'copy',
            pattern: "*.bai"
        ],
        // Gene counts
        [
            path: "${params.proj_dir()}/04.STAR/",
            mode: 'copy',
            pattern: "*.ReadsPerGene.out.tab"
        ],
        // Logs
        [
            path: "${params.proj_dir()}/07.Logs/",
            mode: 'copy',
            pattern: "*.error.log"
        ]
    ]
}
```

### 13.4 Mixing and Collecting Multiple Outputs

**Real Example from main.nf:**
```groovy
// Collect all QC outputs for MultiQC
multiqc_ch = Channel.empty()
    // FastQC outputs
    .mix(FASTQC_RAW.out.fastqc_zip)
    
    // SALMON outputs (extract directory from tuple)
    .mix(SALMON_QUANT.out.salmon_quant.map { it[1] })
    
    // STAR outputs
    .mix(STAR_ALIGN.out.gene_counts)
    .mix(STAR_ALIGN.out.sj_tab)
    .mix(STAR_ALIGN.out.star_log)
    
    // RSeQC outputs
    .mix(RSEQC.out.rseqc_logs)
    
    // Wait for ALL samples to complete
    .collect()

MULTIQC(multiqc_ch)
```

### 13.5 File Pairing and Grouping

**Real Example from validate.nf:**
```groovy
// Map R1 files to create [sample_id, [R1, R2]] tuples
R1_FASTQS_ch = Channel.fromList(all_r1_files)

grouped_samples_ch = R1_FASTQS_ch.map { r1 ->
    // Extract sample ID
    def idx = r1.name.lastIndexOf(Read1_TAG)
    def sample_id = (idx != -1) ? r1.name.take(idx) : r1.simpleName
    
    if (MODE == "PAIRED_END") {
        // Build R2 filename using reverse trick
        def r2_name = r1.name.reverse()
                             .replaceFirst(Read1_TAG.reverse(), Read2_TAG.reverse())
                             .reverse()
        
        // Find R2 in same directory
        def r2 = r1.parent.resolve(r2_name)
        
        if (!r2.exists()) {
            error "R2 pair not found: ${r2_name}"
        }
        
        return [sample_id, [r1, r2]]
    } else {
        return [sample_id, [r1]]
    }
}
```

---

## 14. Best Practices & Tips

### 14.1 Channel vs Groovy Collections

```groovy
// ❌ WRONG: Groovy .collect() on channel
ch = Channel.fromPath("*.fastq")
    .collect { it.name }  // ERROR: Not Groovy method!

// ✅ CORRECT: Use channel operators
ch = Channel.fromPath("*.fastq")
    .map { it.name }

// ❌ WRONG: Creating channel in config
// nextflow.config
params.files = Channel.fromPath("*.fastq")  // ERROR!

// ✅ CORRECT: Create channel in workflow
// main.nf
workflow {
    files_ch = Channel.fromPath(params.pattern)
}
```

### 14.2 Closures for Dynamic Paths

```groovy
// ❌ WRONG: Immediate evaluation
params.project = "test"
params.proj_dir = "${params.base_dir}/${params.project}"
// proj_dir is fixed as ".../test"

// Profile changes params.project to "Xinyi"
// But proj_dir still equals ".../test" ❌

// ✅ CORRECT: Delayed evaluation
params.proj_dir = { "${params.base_dir}/${params.project}" }
// Must call with parentheses: params.proj_dir()
```

### 14.3 Error Handling Best Practices

```groovy
// ✅ GOOD: Fail fast with clear messages
if (!file(ref_genome).exists()) {
    error """
    ❌ Reference genome not found!
    Expected: ${ref_genome}
    Please check your configuration.
    """
}

// ✅ GOOD: Comprehensive error logging
script:
def LOG = "${sample_id}.PROCESS.error.log"
"""
tool --input ${input} \
    1>> "${LOG}" 2>&1 \
    || { echo "❌ ERROR: Tool failed" | tee -a "${LOG}" >&2; exit 1; }

echo "✅ SUCCESS: Tool completed" >> "${LOG}"
"""
```

### 14.4 Resource Allocation Strategy

```groovy
// ✅ GOOD: Use labels for tiers
process {
    withLabel: 'process_low' {
        cpus   = 2
        memory = { 3.GB * task.attempt }
    }
    
    withLabel: 'process_high' {
        cpus   = 8
        memory = { 48.GB * task.attempt }
    }
}

// ✅ GOOD: Memory scales with retries
// Attempt 1: 48GB
// Attempt 2: 96GB (automatic retry if OOM)
```

### 14.5 Profile Organization

```groovy
profiles {
    // Project profiles
    xinyi {
        params.project = "Xinyi"
        params.species = "Human"
    }
    
    // Executor profiles
    sge {
        process.executor = 'sge'
        workDir = "/scratch/work"
    }
    
    // Combine: nextflow run -profile xinyi,sge
}
```

---

## 15. Debugging & Troubleshooting

### 15.1 Debug Channels with .view()
```groovy
workflow {
    fastq_ch = Channel.fromPath("*.fastq")
        .view()  // Print each file
        .map { [it.simpleName, it] }
        .view { "Mapped: $it" }  // Print tuples
}

// Output:
// /path/to/sample1.fastq
// Mapped: [sample1, /path/to/sample1.fastq]
```

### 15.2 Inspect Work Directory
```bash
# Find work directory for specific sample
nextflow log

# Navigate to work directory
cd work/3a/517fad7a2216eb95706e660bc3cc68/

# Check files
ls -lh  # See symlinks
cat .command.sh  # See executed script
cat .command.log  # See stdout/stderr
cat .command.err  # See errors
cat .exitcode  # See exit code
```

### 15.3 Common Errors

#### Error: "No such variable: it"
```groovy
// ❌ WRONG: No implicit 'it' in map with explicit param
.map { file -> it.name }  // ERROR!

// ✅ CORRECT: Use named parameter
.map { file -> file.name }

// ✅ CORRECT: Or use implicit 'it'
.map { it.name }
```

#### Error: "Process terminated for an unknown reason"
```bash
# Check:
# 1. Memory limit exceeded (exit code 137)
# 2. Time limit exceeded
# 3. Disk quota exceeded

# View exit code:
cat work/xx/yyyy.../.exitcode
```

#### Error: "Missing output file(s)"
```groovy
// ❌ WRONG: Pattern doesn't match actual output
output:
path "*.bam"

script:
"""
touch output.BAM  # Different extension!
"""

// ✅ CORRECT: Match exact pattern
output:
path "*.{bam,BAM}"  # Or ensure lowercase in script
```

### 15.4 Resume from Failure
```bash
# First run fails
nextflow run main.nf -profile xinyi

# Fix issue, then resume (uses cache)
nextflow run main.nf -profile xinyi -resume

# Force fresh run (ignore cache)
nextflow run main.nf -profile xinyi

# Clean work directory
nextflow clean -f
```

---

## 16. Quick Reference Table

### 16.1 Channel Factories

| Command | Purpose | Example |
|---------|---------|---------|
| `Channel.fromPath()` | Files matching pattern | `Channel.fromPath("*.fastq")` |
| `Channel.fromFilePairs()` | Auto-pair R1/R2 | `Channel.fromFilePairs("*_{1,2}.fq")` |
| `Channel.value()` | Constant/singleton | `Channel.value(file("ref.fa"))` |
| `Channel.of()` | Literal values | `Channel.of(1, 2, 3)` |
| `Channel.fromList()` | From Groovy list | `Channel.fromList(my_list)` |
| `Channel.empty()` | Empty channel | `Channel.empty()` |

### 16.2 Channel Operators

| Operator | Purpose | Input → Output |
|----------|---------|----------------|
| `.map()` | Transform items | `[1,2,3]` → `[2,4,6]` |
| `.filter()` | Keep matches | `[1,2,3,4]` → `[2,4]` |
| `.collect()` | Wait for all | `1,2,3` → `[1,2,3]` |
| `.flatten()` | Expand lists | `[1,2],[3,4]` → `1,2,3,4` |
| `.unique()` | Remove dupes | `1,2,2,3` → `1,2,3` |
| `.join()` | Merge by key | `[a,1],[a,2]` → `[a,1,2]` |
| `.mix()` | Combine channels | `ch1+ch2` → mixed |
| `.view()` | Debug print | Shows items |

### 16.3 Process Directives

| Directive | Purpose | Example |
|-----------|---------|---------|
| `tag` | Display name | `tag "$sample_id"` |
| `label` | Resource group | `label 'process_high'` |
| `cpus` | CPU cores | `cpus 8` |
| `memory` | RAM allocation | `memory '32 GB'` |
| `time` | Max runtime | `time '2h'` |
| `container` | Docker/Singularity | `container 'ubuntu:20.04'` |
| `publishDir` | Output location | `publishDir '/results'` |
| `errorStrategy` | Error handling | `errorStrategy 'retry'` |
| `maxRetries` | Retry attempts | `maxRetries 3` |
| `when` | Conditional run | `when: params.run == true` |

### 16.4 Input/Output Qualifiers

| Type | Purpose | Example |
|------|---------|---------|
| `val` | Simple value | `val sample_id` |
| `path` | File/directory | `path fastq` |
| `tuple` | Grouped values | `tuple val(id), path(file)` |
| `env` | Environment var | `env MODE` |
| `emit` | Named output | `emit: result` |
| `optional` | May not exist | `optional: true` |

---

## 17. Complete Workflow Example

**Directory Structure:**
```bash
rnaseq_pipeline/
├── main.nf
├── nextflow.config
└── modules/
    ├── validate.nf
    ├── fastqc.nf
    ├── star_align.nf
    └── multiqc.nf
```

**main.nf:**
```groovy
#!/usr/bin/env nextflow
nextflow.enable.dsl=2

include { VALIDATE_INPUT } from './modules/validate.nf'
include { FASTQC } from './modules/fastqc.nf'
include { STAR_ALIGN } from './modules/star_align.nf'
include { MULTIQC } from './modules/multiqc.nf'

workflow {
    // Validate input
    VALIDATE_INPUT(params.raw_fastq_dir())
    sample_ch = VALIDATE_INPUT.out.samples
    
    // QC
    FASTQC(sample_ch)
    
    // Align
    ref_ch = Channel.value(file(params.star_index))
    STAR_ALIGN(sample_ch, ref_ch)
    
    // Aggregate reports
    qc_ch = Channel.empty()
        .mix(FASTQC.out.zip)
        .mix(STAR_ALIGN.out.log)
        .collect()
    
    MULTIQC(qc_ch)
}
```

**nextflow.config:**
```groovy
params {
    project = "test"
    base_dir = "$HOME/scratch"
    proj_dir = { "${params.base_dir}/${params.project}" }
    raw_fastq_dir = { "${params.proj_dir()}/01.FastQ/raw" }
    star_index = "/ref/star_index"
}

profiles {
    standard {
        process.executor = 'local'
    }
    
    sge {
        process.executor = 'sge'
        process.queue = 'all.q'
        workDir = "/scratch/work"
    }
}

process {
    errorStrategy = 'retry'
    maxRetries = 3
    
    withLabel: 'process_low' {
        cpus = 2
        memory = { 3.GB * task.attempt }
    }
}
```

**Run:**
```bash
nextflow run main.nf -profile sge -resume
```

---

## 18. Summary: Nextflow vs Groovy

| Feature | Nextflow | Groovy |
|---------|----------|--------|
| **Purpose** | Workflow orchestration | Logic & validation |
| **When** | Runtime (process execution) | Before runtime (setup) |
| **Data** | Channels (asynchronous pipes) | Lists/Maps (in-memory) |
| **Files** | Channel.fromPath() → stream | files() → immediate Set |
| **Transform** | .map{} on channels | .collect{} on lists |
| **Filter** | .filter{} on channels | .findAll{} on lists |
| **Collect** | Wait for all → emit list | Transform each item |
| **Variables** | params.x (global) | def x (local) |
| **Location** | main.nf, workflows | nextflow.config, validation |

---

## 🎓 Learning Path

1. **Start**: Understand channels vs Groovy collections
2. **Next**: Master process inputs/outputs
3. **Then**: Learn publishDir and configuration
4. **Advanced**: Implement error handling and retry logic
5. **Pro**: Optimize with executors and containers

**Remember**: Nextflow handles **when** and **where** processes run. Groovy handles **what** data to process and **how** to validate it.