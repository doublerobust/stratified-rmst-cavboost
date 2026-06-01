# MoE-K Simulation Environment
# Build:   docker build -t moe-k-sim .
# Run:     docker run --rm -v /path/to/results:/app/results moe-k-sim
# On Omen with GPU:   docker run --rm -v C:\Users\...\results:/app/results moe-k-sim

FROM rocker/r-base:4.6.0

RUN apt-get update && apt-get install -y --no-install-recommends \
    python3 python3-pip python3-venv \
    libcurl4-openssl-dev libssl-dev libxml2-dev \
    git \
    && rm -rf /var/lib/apt/lists/*

# Install R packages
RUN R -e 'install.packages(c( \
    "mvtnorm", "glmnet", "survival", "xgboost", "ranger", \
    "dplyr", "moments", "data.table", "ggplot2", "gridExtra", \
    "RColorBrewer", "viridisLite", "devtools", "furrr", "future", \
    "Matrix" \
  ), repos = "https://cloud.r-project.org", Ncpus = 4)'

# Install Python packages for parallel launcher scripts
RUN python3 -m pip install --no-cache-dir numpy pandas

WORKDIR /app

# Default: just check environment
CMD ["R", "-e", "cat(sessionInfo()$R.version$version.string, '\nPackages:\n'); for(p in c('glmnet','ranger','survival','xgboost','mvtnorm','dplyr','data.table','furrr')) cat('  ', p, ':', as.character(packageVersion(p)), '\n', sep='')"]
