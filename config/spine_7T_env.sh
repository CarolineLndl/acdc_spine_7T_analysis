#!/bin/bash
# =============================================================================
# Environment setup script for spine_7T analysis
# =============================================================================
# This script sets up the required environment for running the spine_7T pipeline.
# It configures temporary directories, PATHs, Python environment, matlab and neuroimaging
# toolboxes (SCT,FSL).
#
# Dependencies:
# - Anaconda/Conda (Python 10 environment: spine_7T_env_py10)
# - Spinal Cord Toolbox (SCT)
# - FSL (FMRIB Software Library)
# - Matlab (for denoising step)
# Note: Adjust paths below according to your filesystem structure.
# =============================================================================

# ----------------------------
# Paths to tools
# ----------------------------
anaconda_dir=$(conda info --base)

SCT_DIR=/cerebro/cerebro1/dataset/spine_7T/derivatives/toolboxes/spinalcordtoolbox
FSLDIR=/cerebro/cerebro1/dataset/bmpd/derivatives/thibault_test/code/toolbox/fsl
MATLAB_DIR=/export01/local/matlab23b

export TMPDIR='/export02/data/tmp'  # temporary files
export MPLBACKEND=Agg              # matplotlib backend for non-GUI rendering

# ----------------------------
# Update PATH
# ----------------------------
export PATH="$SCT_DIR/bin:$PATH"   # spinalcordtoolbox
export PATH=${FSLDIR}/bin:${PATH} # FSL


# ----------------------------
# Python setup
# ----------------------------
source ${anaconda_dir}/etc/profile.d/conda.sh
conda activate spine_7T_env_py10
echo "++ Python executable: $(which python)"
echo "++ Python version: $(python --version)"


# ----------------------------
# FSL setup
# ----------------------------
export PATH=${FSLDIR}/bin:${PATH} # FSL
export FSLDIR PATH
. ${FSLDIR}/etc/fslconf/fsl.sh


# ----------------------------
# Matlab setup (for denoising step)
# ----------------------------
# Uncomment only when running denoising analyses with MATLAB.
# Leaving this commented prevents conflicts between MATLAB’s Qt libraries and SCT’s Qt version.
#LD_PREFIX="${MATLAB_DIR}/sys/os/glnxa64:/cerebro/cerebro1/dataset/bmpd/derivatives/thibault_test/code/toolbox/libraries"
#export  LD_LIBRARY_PATH=/export01/local/matlab23b/bin/glnxa64/



