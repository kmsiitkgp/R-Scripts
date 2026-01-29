#!/bin/bash

# =========================================================================================
# NEXTFLOW PIPELINE LAUNCHER
# =========================================================================================
# This script:
#   1. Reads paths from project_info.yaml
#   2. Auto-detects physical paths (resolving symlinks)
#   3. Builds Singularity bind mounts dynamically
#   4. Launches the Nextflow pipeline with proper configuration
#
# Usage:
#   bash run_nextflow.sh
#
# Prerequisites:
#   - project_info.yaml configured with correct paths
#   - Nextflow and Singularity modules available
#   - Access to HPC cluster (if using -profile sge)
#
# Resume a failed run:
#   Just run the script again - Nextflow automatically resumes from last checkpoint
# =========================================================================================

# -----------------------------------------------------------------------------------------
# CONFIGURATION
# -----------------------------------------------------------------------------------------

# Pipeline directory (where main.nf, project_info.yaml, modules/ are located)
NF_DIR="$HOME/projects/nf_rnaseq"
YAML="${NF_DIR}/project_info.yaml"

# -----------------------------------------------------------------------------------------
# YAML PREPROCESSING
# -----------------------------------------------------------------------------------------
# YAML parsers don't handle tabs well - convert to spaces
# This prevents "mapping values are not allowed here" errors

# Replace tabs with 4 spaces
sed -i 's/\t/    /g' "${YAML}"

# Verify no tabs remain (should print nothing if clean)
grep $'\t' "${YAML}"

# -----------------------------------------------------------------------------------------
# EXTRACT PATHS FROM YAML
# -----------------------------------------------------------------------------------------
# Parse YAML to get directory paths
# These may be logical paths (symlinks) that need to be resolved for Singularity

REF_DIR=$(grep "^[[:space:]]*ref_dir[[:space:]]*:" $YAML | cut -d':' -f2- | xargs | tr -d '"')
READ_DIR=$(grep "^[[:space:]]*read_dir[[:space:]]*:" $YAML | cut -d':' -f2- | xargs | tr -d '"')
BASE_DIR=$(grep "^[[:space:]]*base_dir[[:space:]]*:" $YAML | cut -d':' -f2- | xargs | tr -d '"')
WORK_DIR=$(grep "^[[:space:]]*work_dir[[:space:]]*:" $YAML | cut -d':' -f2- | xargs | tr -d '"')
CACHE_DIR=$(grep "^[[:space:]]*cache_dir[[:space:]]*:" $YAML | cut -d':' -f2- | xargs | tr -d '"')

# -----------------------------------------------------------------------------------------
# RESOLVE PHYSICAL PATHS
# -----------------------------------------------------------------------------------------
# Singularity needs physical paths (not symlinks) for bind mounts
# Example: /scratch -> /gpfs/scratch01 (actual physical location)

P_REF=$(readlink -f "$REF_DIR")
P_READ=$(readlink -f "$READ_DIR")
P_BASE=$(readlink -f "$BASE_DIR")
P_WORK=$(readlink -f "$WORK_DIR")
P_CACHE=$(readlink -f "$CACHE_DIR")

# -----------------------------------------------------------------------------------------
# BUILD SINGULARITY BIND MOUNTS
# -----------------------------------------------------------------------------------------
# Strategy: Bind root directories (e.g., /hpc, /scratch) instead of individual paths
# This covers all subdirectories and avoids complex bind mount conflicts
#
# Why this matters:
#   - If /scratch is a symlink to /gpfs/scratch01, both need to be accessible
#   - Binding parent directories ensures all child paths work
#   - Prevents "No such file or directory" errors inside containers

# Collect all paths (logical and physical)
PATHS_TO_BIND="$REF_DIR $READ_DIR $BASE_DIR $WORK_DIR $CACHE_DIR $P_REF $P_READ $P_BASE $P_WORK $P_CACHE"

# Extract unique root directories (first 2 levels: /hpc, /scratch, /gpfs, etc.)
UNIQUE_ROOTS=$(for p in $PATHS_TO_BIND; do echo "$p" | cut -d'/' -f1-2; done | sort -u | grep '^/')

# Use root directories as bind paths (simplest and most robust approach)
ALL_BIND_PATHS="$UNIQUE_ROOTS"

# Alternative (if you need more specific binds):
# ALL_BIND_PATHS="$UNIQUE_ROOTS $PATHS_TO_BIND"

# Build Singularity bind flags: --bind /hpc --bind /scratch --bind /gpfs
BIND_FLAGS=""
for path in $ALL_BIND_PATHS; do
    if [ -d "$path" ]; then
        BIND_FLAGS+="--bind $path "
    fi
done

echo "Singularity bind flags: $BIND_FLAGS"

# -----------------------------------------------------------------------------------------
# DISPLAY PATH MAPPINGS
# -----------------------------------------------------------------------------------------
# Show how paths map between host (physical) and container (logical)
# Useful for debugging "file not found" issues

echo "---------------------------------------------------------------------------------------------------------------"
printf "%-16s : %-60s -> %s\n" "Mapping:" "PHYSICAL (Host)" "LOGICAL (Container)"
echo "---------------------------------------------------------------------------------------------------------------"
printf "%-16s : %-60s -> %s\n" "Ref Genomes"    "$P_REF"  "$REF_DIR"
printf "%-16s : %-60s -> %s\n" "Raw Reads"      "$P_READ" "$READ_DIR"
printf "%-16s : %-60s -> %s\n" "Project Base"   "$P_BASE" "$BASE_DIR"
printf "%-16s : %-60s -> %s\n" "Temp Work Dir"  "$P_WORK" "$WORK_DIR"
printf "%-16s : %-60s -> %s\n" "Image Cache"    "$P_CACHE" "$CACHE_DIR"
echo "---------------------------------------------------------------------------------------------------------------"

# -----------------------------------------------------------------------------------------
# CLEANUP NEXTFLOW FILES
# -----------------------------------------------------------------------------------------
# Remove problematic characters that cause parsing errors

# Remove UTF-8 BOM (Byte Order Mark) from all Nextflow files
# BOM appears as: 0xEF 0xBB 0xBF at file start (from some text editors)
# Causes error: "Invalid character at start of file"
sed -i '1s/^\xEF\xBB\xBF//' ${NF_DIR}/*.nf ${NF_DIR}/*.config
sed -i '1s/^\xEF\xBB\xBF//' ${NF_DIR}/modules/*.nf

# Clean up whitespace and tabs in module files
# - Strip trailing whitespace (prevents unnecessary git diffs)
# - Convert tabs to 4 spaces (consistent indentation)
sed -i 's/[[:space:]]*$//; s/\t/    /g' ${NF_DIR}/modules/*.nf

# -----------------------------------------------------------------------------------------
# LOAD REQUIRED MODULES
# -----------------------------------------------------------------------------------------
# Load Nextflow and Singularity from HPC environment modules
# Versions may vary by cluster - adjust as needed

module load nextflow/24.10.5
module load singularity-apptainer/1.1.8

# Optional: Clean previous failed runs (use with caution!)
# nextflow clean -f

# -----------------------------------------------------------------------------------------
# EXECUTE NEXTFLOW PIPELINE
# -----------------------------------------------------------------------------------------

# Move to work directory before running
# This prevents Nextflow from cluttering the current directory with temporary files
# Temporary files created: _nf_config_*, .nextflow.log, .nextflow/
mkdir -p "${WORK_DIR}"
cd "${WORK_DIR}"

# Run the pipeline
# -params-file: Load configuration from YAML
# -profile: Use SGE executor (see nextflow.config for other profiles)
# -resume: Automatically resume from last successful checkpoint
# --dynamic_binds: Custom parameter passing bind flags to Singularity (see nextflow.config)
nextflow run "${NF_DIR}/main.nf" \
  -params-file "${YAML}" \
  -profile sge \
  -resume \
  --dynamic_binds "$BIND_FLAGS"
  #-preview

# -----------------------------------------------------------------------------------------
# NOTES ON SINGULARITY BIND MOUNT CHALLENGES
# -----------------------------------------------------------------------------------------
# Problem: Nextflow needs to pass bind mounts to Singularity, but has limited options
#
# Approaches that DIDN'T work:
#   ✗ export NXF_SINGULARITY_RUNOPS="$BIND_FLAGS"
#   ✗ --singularity.runOptions "${BIND_FLAGS}"
#   ✗ --singularity.runOptions "'${BIND_FLAGS}'"
#   ✗ Hardcoding in nextflow.config: runOptions = "--bind /scratch ..."
#
# Solution that WORKS:
#   ✓ Pass bind flags as custom parameter: --dynamic_binds "$BIND_FLAGS"
#   ✓ In nextflow.config: runOptions = "${params.dynamic_binds}"
#
# Why this works:
#   - Custom parameters are evaluated at runtime
#   - Allows dynamic path detection instead of hardcoding
#   - Script calculates binds, passes to Nextflow, Nextflow passes to Singularity

# -----------------------------------------------------------------------------------------
# DEBUGGING TIPS
# -----------------------------------------------------------------------------------------

# Check for UTF-8 BOM in files:
# head -c 3 "${NF_DIR}/main.nf" | od -c
# If output shows: 0000000 357 273 277 # ! /
#   → BOM present, needs removal
# If output shows: 0000000 # ! /
#   → No BOM, file is clean

# View which work directories correspond to which processes:
# nextflow log <run_name> -f name,workdir,status,duration

# List all runs:
# nextflow log

# Clean up old run metadata:
# nextflow clean -f

# Remove temporary Nextflow config directories (safe after pipeline completes):
# rm -rf ~/_nf_config_*

# -----------------------------------------------------------------------------------------
# COMMON ISSUES AND SOLUTIONS
# -----------------------------------------------------------------------------------------

# Issue: "No such file or directory" inside container
# Solution: Check path mappings printed above, ensure physical paths are bound

# Issue: "mapping values are not allowed here" in YAML
# Solution: Run this script - it fixes tabs in YAML automatically

# Issue: Pipeline won't resume after failure
# Solution: Don't run 'nextflow clean' - it deletes cache needed for -resume

# Issue: Out of disk space
# Solution: Clear old work directories: rm -rf ${WORK_DIR}/*
#           Warning: This prevents -resume for those runs!

# Issue: Singularity can't download images
# Solution: Check internet connection, verify cache_dir is writable
