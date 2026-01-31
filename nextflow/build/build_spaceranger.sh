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
## 4) Move to directory containing spaceranger.tar.gz file and spaceranger.def file 
## cd "/mnt/c/Users/kailasamms/OneDrive - Cedars-Sinai Health System/Desktop/build"
## 5) Verify if you can see the files
## ls
## 6) Build Image
## sudo apptainer build spaceranger-4.0.1.sif spaceranger.def
## 7) Upload to singularity_cache dir

# OPTION 2: BUILDING ON HPC?AWS WITH ROOT ACCESS

#!/bin/bash
set -e  # Exit on error

# 1. Move to script directory
cd "$(dirname "$(readlink -f "$0")")"

# 2. Variables & Logging
TAR_FILE="spaceranger-4.0.1.tar.gz"
DOWNLOAD_URL="https://cf.10xgenomics.com/releases/spatial-exp/spaceranger-4.0.1.tar.gz?Expires=1769847544&Key-Pair-Id=APKAI7S6A5RYOXBWRPDA&Signature=ItDbc~hR6b7zLq6mL3CTIpuAoFSrmd4eOYU-E7p2-z21aVwBb~bhIz5l~xGclh8vzr4WsHD7dXk45hdqbxSFCh1CrI9a~9pjAdfISXAkbwM1fVCpmySWwhBouCpccD~Aw~4lcZZ4YqHz0QhXEBMD8EdOHcFY2yUiXeq7poFjkGmr9ILPGdngnggSQ1G2p7-O9Nk981~uT1PHdoVhMzOdfqdCpCbLEVwbgOfxnxlYVNDdhttbCJ50~VXv~2-KygEzzAKJ8--ziQrmRqdYnDlqEgXAt8e3xdsu~PvKTSVehk5C~Eck2AqkKSxqXjziwzuvFd5daHwOZ5feozFEvodjNg__"
IMG_NAME="spaceranger_4.0.1.sif"
IMG_SCRIPT="spaceranger.def"
LOG_FILE="build_log_$(date +%Y%m%d_%H%M%S).txt"

# Start logging everything to both the screen and a file
exec &> >(tee -a "${LOG_FILE}")

echo "=== Build Started at $(date) ==="

# 3. Environment Prep
module load singularity-apptainer/1.1.8 || module load singularity || module load apptainer
export SINGULARITY_TMPDIR=$PWD

# 4. Download Logic
if [ ! -f "${TAR_FILE}" ]; then
    echo "Downloading Space Ranger..."
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
if singularity exec "${IMG_NAME}" spaceranger --version; then
    echo "CHECK PASSED: Space Ranger is alive and well."
    echo "------------------------------------------------"
    
    # 7. Cleanup
    echo "Cleaning up large tarball: ${TAR_FILE}"
    rm "${TAR_FILE}"
    echo "Build process complete. Image: ${IMG_NAME}"
else
    echo "------------------------------------------------"
    echo "CRITICAL ERROR: Space Ranger failed the version check!"
    echo "Keeping ${TAR_FILE} for debugging."
    echo "Check ${LOG_FILE} for details."
    exit 1
fi

echo "=== Build Finished Successfully at $(date) ==="