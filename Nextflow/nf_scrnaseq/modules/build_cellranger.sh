# OPTION 1: BUILDING ON WINDOWS DESKTOP/LAPTOP

## 1) Windows -> Powershell -> wsl --install -d Ubuntu
## 2) Reboot computer
## 3) Update wsl
## Windows -> wsl -> 
## 		sudo apt update && sudo apt upgrade -y
##		sudo apt install -y software-properties-common
##		sudo add-apt-repository -y ppa:apptainer/ppa
##		sudo apt update
##		sudo apt install -y apptainer
## 4) Move to directory containing cellranger.tar.gz file and cellranger.def file 
## cd "/mnt/c/Users/kailasamms/OneDrive - Cedars-Sinai Health System/Desktop/nf_scrnaseq"
## 5) Verify if you can see the files
## ls
## 6) Build Image
## sudo apptainer build cellranger.sif cellranger.def
## 7) Upload to singularity_cache dir

# OPTION 2: BUILDING ON HPC?AWS WITH ROOT ACCESS

#!/bin/bash
set -e  # Exit on error

# 1. Move to script directory
cd "$(dirname "$(readlink -f "$0")")"

# 2. Variables & Logging
TAR_FILE="cellranger-10.0.0.tar.gz"
DOWNLOAD_URL="https://cf.10xgenomics.com/releases/cell-exp/cellranger-10.0.0.tar.gz?Expires=1769765821&Key-Pair-Id=APKAI7S6A5RYOXBWRPDA&Signature=YerBr4hrgtShcJX-IejhWom7wwOStFD-GjakodZ7UlrMLAS9mfQeLmfEx3TNbQG2op7DSJnZW3qbUUxQdMu7apMyeQz5iK36ZsJcrk5TAeMZX8y2k9oxZyIB6xowDjzMMIICd14oOZJcQNpsSuClbvrPyThGQW3JHa4MBrJJKp~-r1Kuq07nOK9ggSg6R9sNd04RLBlpBAjd2Y6~kxcGnbsuD8jIAmBPGpGSkVZbhsjIIfRWdTSEs-nObwfU6w-0BXwsTq~0~saVfVTLIfB03lg3UXu7sE3FHGJ4gn7rN9lfrjqsl6ezBuIJix3lPa34lUX8CAC2GSCSTBj8S3H2Lg__"
IMG_NAME="cellranger_10.0.0.sif"
IMG_SCRIPT="cellranger.def"
LOG_FILE="build_log_$(date +%Y%m%d_%H%M%S).txt"

# Start logging everything to both the screen and a file
exec &> >(tee -a "${LOG_FILE}")

echo "=== Build Started at $(date) ==="

# 3. Environment Prep
module load singularity-apptainer/1.1.8 || module load singularity || module load apptainer
export SINGULARITY_TMPDIR=$PWD

# 4. Download Logic
if [ ! -f "${TAR_FILE}" ]; then
    echo "Downloading Cell Ranger..."
    wget -c -O "${TAR_FILE}" "${DOWNLOAD_URL}"
else
    echo "Using existing tarball: ${TAR_FILE}"
fi

# 5. Build Image
echo "Starting build process..."
singularity build --fakeroot "${IMG_NAME}" "${IMG_SCRIPT}"

# 6. VALIDATION GATE: The "Before Cleaning" Check
echo "------------------------------------------------"
echo "VERIFYING IMAGE INTEGRITY..."
if singularity exec "${IMG_NAME}" cellranger count --version; then
    echo "CHECK PASSED: Cell Ranger is alive and well."
    echo "------------------------------------------------"
    
    # 7. Cleanup
    echo "Cleaning up large tarball: ${TAR_FILE}"
    rm "${TAR_FILE}"
    echo "Build process complete. Image: ${IMG_NAME}"
else
    echo "------------------------------------------------"
    echo "CRITICAL ERROR: Cell Ranger failed the version check!"
    echo "Keeping ${TAR_FILE} for debugging."
    echo "Check ${LOG_FILE} for details."
    exit 1
fi

echo "=== Build Finished Successfully at $(date) ==="