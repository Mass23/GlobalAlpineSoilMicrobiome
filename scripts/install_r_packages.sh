#!/usr/bin/env bash
# Install CRAN and GitLab R packages not available on conda-forge into the activated conda env.
# Usage (from Snakemake this runs inside the created conda env):
#   conda activate globalalpine-r-download
#   bash scripts/install_r_packages.sh

set -euo pipefail

# Ensure remotes is available
Rscript -e "if(!'remotes' %in% installed.packages()[,'Package']) install.packages('remotes', repos='https://cloud.r-project.org')"

# Install CopernicusDEM from CRAN (required for DEM helpers)
Rscript -e "install.packages('CopernicusDEM', repos='https://cloud.r-project.org')"

# Install rchelsa from GitLab (karger/rchelsa)
Rscript -e "remotes::install_gitlab('karger/rchelsa')"

cat <<'MSG'
Installed CRAN packages: CopernicusDEM and GitLab package rchelsa
If you need additional packages, add them to this script.
MSG
