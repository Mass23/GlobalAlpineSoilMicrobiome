#!/usr/bin/env bash
# Install CRAN R packages not available on conda-forge into the activated conda env.
# Usage:
#   conda activate globalalpine-r
#   bash scripts/install_r_packages.sh

set -euo pipefail

Rscript -e "install.packages(c('CopernicusDEM'), repos='https://cloud.r-project.org')"

cat <<'MSG'
Installed CRAN packages: CopernicusDEM
If you need other CRAN packages, add them to the Rscript call above.
MSG
