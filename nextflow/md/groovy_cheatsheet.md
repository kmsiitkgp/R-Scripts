# 🎯 Complete Groovy Cheatsheet for Nextflow Users

> **Purpose**: Groovy handles all logic, validation, and data manipulation in Nextflow pipelines  
> **Key Concept**: Groovy runs BEFORE workflow execution for setup and validation

---

## 📋 Table of Contents
1. [Variables & Scoping](#1-variables--scoping)
2. [Collections (Lists, Sets, Maps)](#2-collections-lists-sets-maps)
3. [Strings & String Interpolation](#3-strings--string-interpolation)
4. [Regular Expressions (Regex)](#4-regular-expressions-regex)
5. [File & Path Operations](#5-file--path-operations)
6. [Closures (Lambda Functions)](#6-closures-lambda-functions)
7. [Control Flow & Logic](#7-control-flow--logic)
8. [Loops & Iteration](#8-loops--iteration)
9. [Data Transformation](#9-data-transformation)
10. [Operators & Comparisons](#10-operators--comparisons)
11. [Type Conversion](#11-type-conversion)
12. [Debugging & Logging](#12-debugging--logging)

---

## 1. Variables & Scoping

### 1.1 Local Variables (def)
```groovy
// LOCAL SCOPE - Only visible in current block
def sample_id = "Sample1"
def cpus = 8
def is_paired = true

// Use 'def' inside workflows, functions, closures
workflow {
    def fastq_files = files("*.fq.gz")  // Only visible in this workflow
}
```

### 1.2 Global Variables (no def)
```groovy
// GLOBAL SCOPE - Visible everywhere in script (AVOID IN DSL2)
sample_id = "Sample1"  // Dangerous! Can be overwritten anywhere

// Only use for script-level constants
VALID_PATTERN = ~/.*((_R|_r)[12]).*\.f(q|astq)\.gz/
```

### 1.3 Pipeline Parameters (params)
```groovy
// BEST PRACTICE - Use params for pipeline-wide configuration
params.project = "Xinyi"
params.species = "Human"
params.cpus = 8

// Access anywhere in pipeline
println "Running project: ${params.project}"

// Override from command line
// nextflow run main.nf --project "NewProject" --cpus 16
```

### 1.4 Final Variables (Constants)
```groovy
// Cannot be reassigned
final String GENOME_VERSION = "GRCh38"
final int MAX_RETRIES = 3

// GENOME_VERSION = "GRCm39"  // ERROR: Cannot reassign
```

### 1.5 Closures for Dynamic Evaluation
```groovy
// WITHOUT closure - evaluates immediately
params.project = "test"
params.proj_dir = "${params.base_dir}/${params.project}"  
// proj_dir is now ".../test" (fixed!)

// WITH closure - evaluates on demand
params.proj_dir = { "${params.base_dir}/${params.project}" }
// Later: profile sets params.project = "Xinyi"
println params.proj_dir()  // NOW evaluates to ".../Xinyi" ✅
```

**Real Example from nextflow.config:**
```groovy
params {
    base_dir = "/scratch"
    project = "test"
    
    // Closure delays evaluation until runtime
    proj_dir = { "${params.base_dir}/${params.expt}/${params.project}" }
    fastq_dir = { "${params.proj_dir()}/01.FastQ" }
    raw_fastq_dir = { "${params.fastq_dir()}/raw" }
}

// Must call with parentheses
println params.proj_dir()      // Returns: "/scratch/RNASeq/Xinyi"
println params.raw_fastq_dir() // Returns: "/scratch/RNASeq/Xinyi/01.FastQ/raw"
```

---

## 2. Collections (Lists, Sets, Maps)

### 2.1 Lists (Ordered, Allows Duplicates)
```groovy
// Creation
def samples = ["Sample1", "Sample2", "Sample3"]
def numbers = [1, 2, 3, 4, 5]
def mixed = ["text", 42, true, 3.14]

// Access by index
println samples[0]        // "Sample1"
println samples[-1]       // "Sample3" (last element)

// Common operations
samples.size()            // 3
samples.isEmpty()         // false
samples.contains("Sample1")  // true
samples.add("Sample4")    // [Sample1, Sample2, Sample3, Sample4]
samples << "Sample5"      // Append operator (same as .add())

// Slicing
samples[0..1]             // [Sample1, Sample2]
samples[1..-1]            // [Sample2, Sample3]

// Sorting
numbers.sort()            // [1, 2, 3, 4, 5]
samples.sort()            // Alphabetical

// Remove duplicates
def dups = [1, 2, 2, 3, 3, 3]
dups.unique()             // [1, 2, 3]

// Reverse
samples.reverse()         // [Sample3, Sample2, Sample1]

// Min/Max
numbers.min()             // 1
numbers.max()             // 5

// Join to string
samples.join(", ")        // "Sample1, Sample2, Sample3"
```

**Real Example from validate.nf:**
```groovy
// Get all FASTQ files and sort them
fastq_files_set = files("${fastq_dir}/*.f*q.gz")
fastq_files_list = fastq_files_set.collect()
fastq_files = fastq_files_list.sort()  // Sorted list

// Count files by type
def total_files = fastq_files.size()
def r1_files = fastq_files.findAll { it.name.contains("_R1") }
def r2_files = fastq_files.findAll { it.name.contains("_R2") }

println "Total: ${total_files}, R1: ${r1_files.size()}, R2: ${r2_files.size()}"
```

### 2.2 Sets (Unordered, No Duplicates)
```groovy
// Creation (converts list to set)
def numbers_set = [1, 2, 2, 3, 3] as Set  // {1, 2, 3}

// From files() - returns a Set to avoid duplicates
def fastq_set = files("*.fq.gz")          // PathSet

// Cannot access by index
// fastq_set[0]  // ERROR!

// Check membership (FAST)
fastq_set.contains(file("sample1.fq.gz"))  // true

// Convert to list for indexing
def fastq_list = fastq_set.collect()
println fastq_list[0]  // Now works
```

### 2.3 Maps (Key-Value Pairs)
```groovy
// Creation
def sample_info = [
    id: "Sample1",
    species: "Human",
    reads: 50000000,
    paired: true
]

// Access
sample_info.id           // "Sample1"
sample_info['species']   // "Human" (alternative syntax)

// Add/Update
sample_info.genome = "GRCh38"
sample_info['version'] = "v1"

// Check key existence
sample_info.containsKey('id')  // true

// Iterate over keys/values
sample_info.each { key, value ->
    println "Key: ${key}, Value: ${value}"
}
```

**Real Example from nextflow.config:**
```groovy
// Reference genome maps
params.ref_fasta = [
    Human : "${params.ref_dir}/Human/Homo_sapiens.GRCh38.fa",
    Mouse : "${params.ref_dir}/Mouse/Mus_musculus.GRCm39.fa"
]

params.ref_gtf = [
    Human : "${params.ref_dir}/Human/Homo_sapiens.GRCh38.113.gtf",
    Mouse : "${params.ref_dir}/Mouse/Mus_musculus.GRCm39.113.gtf"
]

// Access by species
def human_fasta = params.ref_fasta[params.species]  // When species="Human"
```

---

## 3. Strings & String Interpolation

### 3.1 String Literals
```groovy
// Single quotes - no interpolation
def name1 = 'Sample1'
def path1 = '/home/user/${name1}'  // Literal: /home/user/${name1}

// Double quotes - with interpolation
def name2 = "Sample2"
def path2 = "/home/user/${name2}"  // Result: /home/user/Sample2

// Triple quotes - multi-line strings
def message = """
This is a multi-line
string with ${name2}
"""

// GStrings (Groovy Strings with interpolation)
def sample_id = "Sample1"
def filename = "${sample_id}_R1.fastq.gz"  // GString
```

### 3.2 String Interpolation
```groovy
def sample = "Sample1"
def replicate = 2
def species = "Human"

// Basic interpolation
def msg1 = "Processing ${sample}"                    // "Processing Sample1"

// Expressions in interpolation
def msg2 = "Rep ${replicate + 1}"                    // "Rep 3"
def msg3 = "Species: ${species.toUpperCase()}"       // "Species: HUMAN"

// Nested property access
def msg4 = "Length: ${sample.length()}"              // "Length: 7"

// Conditional in interpolation
def msg5 = "Type: ${replicate > 1 ? 'Multi' : 'Single'}"  // "Type: Multi"
```

**Real Example from rseqc.nf:**
```groovy
// Dynamic log filename
def LOG = "${sample_id}.RSEQC.error.log"

// Command with interpolation
"""
echo "✅ SUCCESS: Gene body coverage calculation completed for ${sample_id}" >> "${LOG}"
"""
```

### 3.3 String Methods
```groovy
def str = "Sample_R1_001.fastq.gz"

// Case conversion
str.toUpperCase()         // "SAMPLE_R1_001.FASTQ.GZ"
str.toLowerCase()         // "sample_r1_001.fastq.gz"

// Trimming
"  Sample1  ".trim()      // "Sample1"

// Splitting
str.split("_")            // ["Sample", "R1", "001.fastq.gz"]
str.tokenize("_")         // Same as split (Groovy method)

// Substring
str.substring(0, 6)       // "Sample"
str.take(6)               // "Sample" (Groovy way)
str.drop(7)               // "R1_001.fastq.gz" (remove first 7 chars)

// Replace
str.replace("R1", "R2")   // "Sample_R2_001.fastq.gz"
str.replaceAll("\\d+", "X")  // "Sample_RX_X.fastq.gz" (regex)

// Contains
str.contains("R1")        // true
str.startsWith("Sample")  // true
str.endsWith(".gz")       // true

// Reverse
"ABC".reverse()           // "CBA"

// Repeat
"A" * 3                   // "AAA"

// Find index
str.indexOf("R1")         // 7
str.lastIndexOf("_")      // 10 (last occurrence)
```

**Real Example from validate.nf:**
```groovy
// Extract sample ID by removing read tag
def r1_name = "Sample1_Tumor_R1.fastq.gz"
def Read1_TAG = "_R1"

def idx = r1_name.lastIndexOf(Read1_TAG)  // Find last occurrence
def sample_id = r1_name.take(idx)          // "Sample1_Tumor"
```

---

## 4. Regular Expressions (Regex)

### 4.1 Creating Regex Patterns
```groovy
// Slashy syntax (recommended)
def pattern1 = /sample_R[12]\.fastq\.gz/

// Tilde syntax (creates Pattern object)
def pattern2 = ~/.*((_R|_r)[12]).*\.f(q|astq)\.gz/

// From string (must escape backslashes)
def pattern3 = "sample_R[12]\\.fastq\\.gz"
```

### 4.2 Matching Operators
```groovy
def filename = "Sample1_R1.fastq.gz"

// ==~ : Full string match
filename ==~ /.*R1.*\.fastq\.gz/     // true (entire string matches)
filename ==~ /R1/                     // false (doesn't match full string)

// =~ : Partial match (returns Matcher object)
filename =~ /R1/                      // true (contains "R1")
def matcher = filename =~ /Sample(\d+)/
if (matcher.find()) {
    println matcher.group(1)          // "1" (captured group)
}

// Find operator
filename.find(/R\d/)                  // "R1" (first match)

// Match operator  
filename.matches(/.*R1.*/)            // true (Java method)
```

**Real Example from validate.nf:**
```groovy
// Validation pattern for FASTQ files
def VALID_PATTERN = ~/.*((_Tumor|_Normal))?.*((_R|_r)[12]).*\.f(q|astq)\.gz/

// Classify files
valid_files = fastq_files.findAll { it.name ==~ VALID_PATTERN }
invalid_files = fastq_files.findAll { !(it.name ==~ VALID_PATTERN) }

// Detect read type (case-sensitive)
r1_files = valid_files.findAll { it.name.contains("_r1") }
R1_files = valid_files.findAll { it.name.contains("_R1") }

// Pattern breakdown:
// .* = any characters (start)
// ((_Tumor|_Normal))? = optional Tumor/Normal designation
// .* = any characters (middle)
// ((_R|_r)[12]) = required _R1/_R2 or _r1/_r2
// .* = any characters before extension
// \.f(q|astq)\.gz = .fq.gz or .fastq.gz extension
```

### 4.3 String Replacement with Regex
```groovy
def str = "Sample_R1_001_L001.fastq.gz"

// Replace all digits with X
str.replaceAll(/\d+/, "X")           // "Sample_RX_X_LX.fastq.gz"

// Replace first occurrence
str.replaceFirst(/R1/, "R2")         // "Sample_R2_001_L001.fastq.gz"

// Capture groups
def result = str.replaceAll(/Sample_R(\d)/, 'Read$1')
// "Read1_001_L001.fastq.gz"
```

**Real Example from validate.nf (R1 to R2 conversion):**
```groovy
// Convert R1 filename to R2 filename using reverse trick
def r1_name = "Sample1_R1_data_R1.fq.gz"
def Read1_TAG = "_R1"
def Read2_TAG = "_R2"

// Reverse the string, replace first occurrence, reverse back
// This prevents replacing "R1" in sample name
def r2_name = r1_name.reverse()
                     .replaceFirst(Read1_TAG.reverse(), Read2_TAG.reverse())
                     .reverse()
// Result: "Sample1_R1_data_R2.fq.gz" (only last _R1 changed)
```

---

## 5. File & Path Operations

### 5.1 File Objects
```groovy
// Create file object
def fasta = new File("/path/to/genome.fa")
def gtf = file("/path/to/annotation.gtf")  // Nextflow helper

// File properties
fasta.name              // "genome.fa"
fasta.baseName          // "genome" (without extension)
fasta.extension         // "fa"
fasta.parent            // "/path/to" (directory)
fasta.absolutePath      // Full path
fasta.size()            // File size in bytes

// File checks
fasta.exists()          // true/false
fasta.isFile()          // true (not a directory)
fasta.isDirectory()     // false
fasta.canRead()         // true/false
fasta.canWrite()        // true/false

// Read file content
fasta.text              // Read entire file as string
fasta.bytes             // Read as byte array
```

**Real Example from prep_reference.nf:**
```groovy
// Extract base name from GTF for BED file naming
def gtf_basename = ref_gtf.baseName 
def ref_bed = "${gtf_basename}.bed"

// Example: Homo_sapiens.GRCh38.113.gtf → Homo_sapiens.GRCh38.113.bed
```

### 5.2 Glob Patterns with files()
```groovy
// files() returns a Set of files matching pattern
def all_fastq = files("*.fastq.gz")
def r1_files = files("*_R1.fastq.gz")
def all_in_dir = files("/path/to/data/*.fq")

// Recursive search with **
def nested_files = files("**/*.fastq.gz")  // All FASTQ in subdirectories

// Multiple patterns
def multi = files("*.{fq,fastq}.gz")       // Both .fq.gz and .fastq.gz

// From variable
def dir = "/data"
def pattern = "${dir}/**/*.fq.gz"
def found = files(pattern)
```

**Real Example from validate.nf:**
```groovy
// Load all FASTQ files from raw directory
fastq_dir = data_dir.toString()
fastq_files_set = files("${fastq_dir}/*.f*q.gz")  // .fq.gz or .fastq.gz
fastq_files_list = fastq_files_set.collect()       // Convert Set to List
fastq_files = fastq_files_list.sort()              // Sort alphabetically
```

### 5.3 Path Manipulation
```groovy
def r1 = file("/data/Sample1_R1.fastq.gz")

// Get parent directory
def parent_dir = r1.parent                 // /data

// Build new path in same directory
def r2 = r1.parent.resolve("Sample1_R2.fastq.gz")

// Check existence
if (!r2.exists()) {
    error "R2 file not found: ${r2}"
}

// toString() for full path
println r1.toString()     // /data/Sample1_R1.fastq.gz
```

**Real Example from validate.nf (Pairing R1 and R2):**
```groovy
grouped_samples_ch = R1_FASTQS_ch.map { r1 ->
    // Extract sample ID
    def idx = r1.name.lastIndexOf(Read1_TAG)
    def sample_id = (idx != -1) ? r1.name.take(idx) : r1.simpleName
    
    if (MODE == "PAIRED_END") {
        // Build R2 filename
        def r1_name = r1.name
        def r2_name = r1_name.reverse()
                             .replaceFirst(Read1_TAG.reverse(), Read2_TAG.reverse())
                             .reverse()
        
        // Find R2 in same directory as R1
        def r2 = r1.parent.resolve(r2_name)
        
        // Verify R2 exists
        if (!r2.exists()) {
            error "R2 pair not found: ${r2_name}"
        }
        
        return [sample_id, [r1, r2]]
    } else {
        return [sample_id, [r1]]
    }
}
```

### 5.4 Reading Files
```groovy
def config_file = new File("config.txt")

// Read entire file
def content = config_file.text

// Read line by line
config_file.eachLine { line ->
    println line
}

// Read with line numbers
config_file.eachLine { line, num ->
    println "${num}: ${line}"
}

// Read into list
def lines = config_file.readLines()

// Check if file exists before reading
if (config_file.exists()) {
    def data = config_file.text
}
```

---

## 6. Closures (Lambda Functions)

### 6.1 Basic Closure Syntax
```groovy
// Closure with explicit parameter
def multiply = { x -> x * 2 }
println multiply(5)  // 10

// Closure with implicit 'it' parameter
def square = { it * it }
println square(4)    // 16

// Closure with multiple parameters
def add = { a, b -> a + b }
println add(3, 7)    // 10

// Closure with no parameters
def greet = { -> println "Hello!" }
greet()              // Hello!
```

### 6.2 Closures with Collections
```groovy
def numbers = [1, 2, 3, 4, 5]

// Transform each element
numbers.each { println it * 2 }

// Collect (map) - create new list
def doubled = numbers.collect { it * 2 }  // [2, 4, 6, 8, 10]

// Filter elements
def evens = numbers.findAll { it % 2 == 0 }  // [2, 4]

// Find first match
def firstEven = numbers.find { it % 2 == 0 }  // 2

// Count matches
def evenCount = numbers.count { it % 2 == 0 }  // 2

// Check if any/all match
numbers.any { it > 3 }   // true
numbers.every { it > 0 } // true
```

**Real Example from validate.nf:**
```groovy
// Filter valid files using closure
valid_files = fastq_files.findAll { it.name ==~ VALID_PATTERN }

// Explicit parameter name for clarity
valid_files = fastq_files.findAll { file -> file.name ==~ VALID_PATTERN }

// Classify files by type
fq_gz_files = valid_files.findAll { it.name.endsWith(".fq.gz") }
fastq_gz_files = valid_files.findAll { it.name.endsWith(".fastq.gz") }
```

### 6.3 Closures for Dynamic Evaluation
```groovy
// Closure stored as variable (delays evaluation)
params.base_dir = "/scratch"
params.project = "test"

// Without closure - evaluates NOW
params.proj_dir = "${params.base_dir}/${params.project}"  
// Fixed as "/scratch/test"

// With closure - evaluates LATER
params.proj_dir = { "${params.base_dir}/${params.project}" }
// Changes when project changes

// Later in code
params.project = "Xinyi"
println params.proj_dir()  // "/scratch/Xinyi" ✅
```

---

## 7. Control Flow & Logic

### 7.1 If-Else Statements
```groovy
def mode = "PAIRED_END"

// Basic if
if (mode == "PAIRED_END") {
    println "Processing paired-end data"
}

// If-else
if (mode == "PAIRED_END") {
    println "PE mode"
} else {
    println "SE mode"
}

// If-else if-else
def read_count = 50000000
if (read_count < 10000000) {
    println "Low coverage"
} else if (read_count < 50000000) {
    println "Medium coverage"
} else {
    println "High coverage"
}

// Ternary operator (inline if-else)
def label = (mode == "PAIRED_END") ? "PE" : "SE"
def status = (read_count > 20000000) ? "Pass" : "Fail"
```

**Real Example from rseqc.nf:**
```groovy
// Translate Nextflow mode to RSeQC parameter
def SEQUENCING_MODE = (mode == "PAIRED_END") ? "PE" : "SE"

// Conditional script execution
if (mode == "PAIRED_END") {
    inner_distance.py \
        --input-file "${bam}" \
        --refgene "${ref_bed}" \
        --out-prefix "${sample_id}"
}
```

### 7.2 Switch Statements
```groovy
def species = "Human"

switch(species) {
    case "Human":
        println "Using GRCh38"
        break
    case "Mouse":
        println "Using GRCm39"
        break
    default:
        println "Unknown species"
}

// Switch with patterns
def filename = "Sample_R1.fastq.gz"
switch(filename) {
    case ~/.*R1.*/:
        println "Read 1 file"
        break
    case ~/.*R2.*/:
        println "Read 2 file"
        break
}
```

### 7.3 Null Safety
```groovy
def value = null

// Safe navigation operator (?.)
def length = value?.length()  // null (doesn't crash)

// Elvis operator (?:) - default value
def name = null
def displayName = name ?: "Unknown"  // "Unknown"

// With safe navigation
def sample = null
def id = sample?.id ?: "default"  // "default"
```

---

## 8. Loops & Iteration

### 8.1 For Loops
```groovy
// Range loops
for (i in 1..5) {
    println i  // 1, 2, 3, 4, 5
}

for (i in 1..<5) {
    println i  // 1, 2, 3, 4 (exclusive end)
}

// List iteration
def samples = ["S1", "S2", "S3"]
for (sample in samples) {
    println sample
}

// With index
for (i in 0..<samples.size()) {
    println "${i}: ${samples[i]}"
}
```

### 8.2 Each (Groovy Way)
```groovy
def samples = ["Sample1", "Sample2", "Sample3"]

// each - side effects only
samples.each { sample ->
    println "Processing ${sample}"
}

// eachWithIndex
samples.eachWithIndex { sample, idx ->
    println "${idx}: ${sample}"
}

// each on maps
def info = [name: "Sample1", reads: 1000000]
info.each { key, value ->
    println "${key} = ${value}"
}
```

**Real Example from validate.nf:**
```groovy
// Print invalid filenames
if (invalid_files.size() > 0) {
    println "ERROR: ${invalid_files.size()} INVALID FILE(S) DETECTED"
    
    invalid_files.each { file ->
        println "  - ${file.name}"
    }
    
    error "Please rename files to match pattern"
}
```

### 8.3 While Loops
```groovy
def count = 0
while (count < 5) {
    println count
    count++
}

// While with complex condition
def files = files("*.fastq")
while (!files.isEmpty()) {
    def file = files[0]
    println file.name
    files = files.drop(1)
}
```

---

## 9. Data Transformation

### 9.1 collect() - Transform Elements
```groovy
def numbers = [1, 2, 3, 4, 5]

// Square each number
def squared = numbers.collect { it * it }  // [1, 4, 9, 16, 25]

// Transform to strings
def strings = numbers.collect { "Number: ${it}" }
// ["Number: 1", "Number: 2", ...]

// Extract properties
def files = files("*.fastq")
def names = files.collect { it.name }
def sizes = files.collect { it.size() }
```

**Real Example from validate.nf:**
```groovy
// Convert Set to List using collect()
fastq_files_set = files("${fastq_dir}/*.f*q.gz")
fastq_files_list = fastq_files_set.collect()  // Set → List
fastq_files = fastq_files_list.sort()
```

### 9.2 findAll() - Filter Elements
```groovy
def numbers = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10]

// Get even numbers
def evens = numbers.findAll { it % 2 == 0 }  // [2, 4, 6, 8, 10]

// Get numbers > 5
def large = numbers.findAll { it > 5 }  // [6, 7, 8, 9, 10]

// Filter files
def fastq_files = files("*.gz")
def r1_files = fastq_files.findAll { it.name.contains("_R1") }
def valid = fastq_files.findAll { it.name ==~ /.*R[12]\.fastq\.gz/ }
```

### 9.3 groupBy() - Group Elements
```groovy
def samples = [
    [id: "S1", type: "Tumor"],
    [id: "S2", type: "Normal"],
    [id: "S3", type: "Tumor"]
]

// Group by type
def grouped = samples.groupBy { it.type }
// [Tumor: [[id:S1, type:Tumor], [id:S3, type:Tumor]], 
//  Normal: [[id:S2, type:Normal]]]

// Group files by extension
def files = files("*.*")
def by_ext = files.groupBy { it.extension }
```

### 9.4 join() - Combine to String
```groovy
def samples = ["Sample1", "Sample2", "Sample3"]

samples.join(", ")        // "Sample1, Sample2, Sample3"
samples.join("\n")        // Multi-line string
samples.join(" ")         // Space-separated

// Join file paths
def paths = ["/data", "project", "fastq"]
paths.join("/")           // "/data/project/fastq"
```

**Real Example from star_align.nf:**
```groovy
// Join STAR arguments from params list
STAR \
    --genomeDir "${star_index_dir}" \
    ${params.STAR_ARGS.join(' ')} \
    --runThreadN "${task.cpus}"

// params.STAR_ARGS = ["--runMode", "alignReads", "--twopassMode", "Basic"]
// Result: --runMode alignReads --twopassMode Basic
```

---

## 10. Operators & Comparisons

### 10.1 Comparison Operators
```groovy
// Equality
5 == 5        // true
5 != 3        // true

// Relational
5 > 3         // true
5 < 10        // true
5 >= 5        // true
5 <= 5        // true

// String comparison
"abc" == "abc"      // true
"abc" < "def"       // true (alphabetically)
```

### 10.2 Logical Operators
```groovy
// AND
true && true        // true
true && false       // false

// OR
true || false       // true
false || false      // false

// NOT
!true               // false
!false              // true

// Combined
(age > 18) && (country == "US")
(status == "Tumor") || (status == "Normal")
```

**Real Example from validate.nf:**
```groovy
// Complex validation conditions
if (MODE == "PAIRED_END" && r1_files.size() * 2 == valid_files.size()) {
    println "Valid paired-end dataset"
}

// Negative validation
invalid_files = fastq_files.findAll { !(it.name ==~ VALID_PATTERN) }
```

### 10.3 Membership Operators
```groovy
def samples = ["S1", "S2", "S3"]

// in operator
"S1" in samples       // true
"S4" in samples       // false

// Regex match
"Sample_R1.fastq" ==~ /.*R1.*/    // true (full match)
"Sample_R1.fastq" =~ /R1/          // true (contains)
```

### 10.4 Arithmetic Operators
```groovy
// Basic math
5 + 3         // 8
5 - 3         // 2
5 * 3         // 15
5 / 2         // 2.5
5.intdiv(2)   // 2 (integer division)
5 % 2         // 1 (modulo/remainder)
5 ** 2        // 25 (power)

// Increment/decrement
def x = 5
x++           // x is now 6
++x           // x is now 7
x--           // x is now 6

// Compound assignment
x += 3        // x = x + 3
x *= 2        // x = x * 2
```

**Real Example from validate.nf:**
```groovy
// Calculate sample counts for paired-end
if (MODE == "PAIRED_END") {
    // .intdiv(2) ensures integer division
    TOTAL_SAMPLES = valid_files.size().intdiv(2)
    N_TUMOR_SAMPLES = tumor_files.size().intdiv(2)
    N_NORMAL_SAMPLES = normal_files.size().intdiv(2)
}
```

### 10.5 Spread Operator (*)
```groovy
def list1 = [1, 2, 3]
def list2 = [4, 5, 6]

// Combine lists
def combined = [*list1, *list2]  // [1, 2, 3, 4, 5, 6]

// Call method on each element
def names = ["alice", "bob"]
def upper = names*.toUpperCase()  // ["ALICE", "BOB"]
```

---

## 11. Type Conversion

### 11.1 Explicit Conversion
```groovy
// String to number
def str = "42"
def num = str.toInteger()        // 42
def dbl = str.toDouble()         // 42.0

// Number to string
def number = 42
def text = number.toString()     // "42"

// List to Set (remove duplicates)
def list = [1, 2, 2, 3]
def set = list as Set            // {1, 2, 3}
def set2 = list.toSet()          // Same

// Set to List
def mySet = [1, 2, 3] as Set
def myList = mySet as List

// Array to List
String[] arr = ["a", "b", "c"]
def list2 = arr as List
```

### 11.2 Parsing
```groovy
// Parse integers
Integer.parseInt("42")           // 42

// Parse doubles
Double.parseDouble("3.14")       // 3.14

// Parse boolean
Boolean.parseBoolean("true")     // true
"true".toBoolean()               // true
```

### 11.3 Type Checking
```groovy
def value = "hello"

value instanceof String          // true
value instanceof Number          // false

// Check if number
value.isNumber()                 // false (Groovy extension)
"123".isNumber()                 // true
```

---

## 12. Debugging & Logging

### 12.1 Print Statements
```groovy
// Basic print
println "Hello, World!"

// Print variable
def sample = "Sample1"
println sample
println "Sample: ${sample}"

// Print without newline
print "Loading... "
println "Done!"

// Formatted output
printf("Sample: %s, Reads: %d%n", sample, 1000000)
```

### 12.2 Inspecting Variables
```groovy
def data = [id: "Sample1", reads: 1000000]

// Dump variable content
println data.dump()

// Inspect (similar to dump)
println data.inspect()

// Print type
println data.getClass()

// Print properties
data.properties.each { key, value ->
    println "${key} = ${value}"
}
```

### 12.3 Logging
```groovy
// Using log object (in Nextflow)
log.info "Pipeline started"
log.warn "Low read count detected"
log.error "File not found"
log.debug "Debug information"

// Formatted messages
log.info """
    Project : ${params.project}
    Species : ${params.species}
    Samples : ${sample_count}
"""
```

**Real Example from main.nf:**
```groovy
log.info """
    ===========================================
    RNA-SEQ PIPELINE
    ===========================================
    Project : ${params.project}
    Species : ${params.species}
    Genome  : ${params.genome_version}
    ===========================================
"""
```

### 12.4 Error Handling
```groovy
// Throw error and stop execution
if (!file("ref.fa").exists()) {
    error "Reference genome not found!"
}

// Try-catch (rarely used in Nextflow)
try {
    def result = 10 / 0
} catch (ArithmeticException e) {
    println "Error: ${e.message}"
}

// Assert
assert sample_count > 0 : "No samples found"
```

**Real Example from validate.nf:**
```groovy
// Fail fast with descriptive error
if (valid_files.size() == 0) {
    println "\n" + "!" * 50
    println " ERROR: NO VALID FASTQ FILES FOUND"
    println "!" * 50
    println " Search Path: ${params.raw_fastq_dir()}"
    println "!" * 50 + "\n"
    error "Aborting: Zero files matched the pattern"
}
```

---

## 📚 Summary: When to Use What

| Task | Groovy Method |
|------|---------------|
| **Store pipeline-wide settings** | `params.variable = value` |
| **Local temporary variables** | `def variable = value` |
| **Delay path evaluation** | Use closure: `{ "path" }` |
| **Work with ordered data** | Use Lists `[1, 2, 3]` |
| **Remove duplicates** | Use Sets `.toSet()` or `as Set` |
| **Key-value pairs** | Use Maps `[key: value]` |
| **Variable in strings** | `"text ${var}"` (double quotes) |
| **Match entire string** | `str ==~ /pattern/` |
| **Check substring** | `str =~ /pattern/` or `.contains()` |
| **Get files matching pattern** | `files("*.fastq")` |
| **Read file content** | `file.text` or `file.readLines()` |
| **Transform each element** | `.collect { it * 2 }` |
| **Filter elements** | `.findAll { it > 5 }` |
| **Loop with side effects** | `.each { println it }` |
| **Early error exit** | `error "message"` |
| **Debug output** | `println` or `log.info` |

---

## 🎓 Best Practices

1. **Use `def` for local scope** - Prevents variable collisions
2. **Use closures for dynamic paths** - Allows profile overrides
3. **Prefer `.collect()` over manual loops** - More functional, cleaner
4. **Use `.findAll()` for filtering** - More readable than manual iteration
5. **Always validate file existence** - Fail fast with `.exists()`
6. **Use regex for pattern matching** - More powerful than simple string matching
7. **Convert Sets to Lists when needed** - Sets don't support indexing
8. **Use descriptive variable names** - `sample_id` not `s`
9. **Log important information** - Use `log.info` for pipeline progress
10. **Fail with clear error messages** - Help users fix issues quickly</parameter>