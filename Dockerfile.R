# Install system dependencies
RUN apt-get update && apt-get install -y \
libcurl4-openssl-dev \
libssl-dev \
libxml2-dev \
git \
unzip \
&& rm -rf /var/lib/apt/lists/*

  # Install Bioconductor and CRAN packages
  RUN R -e "install.packages(c('devtools', 'ggplot2', 'stringr', 'igraph', 'Matrix', 'rlang'))"
RUN R -e "install.packages('BiocManager')"
RUN R -e "BiocManager::install(c('DropletUtils', 'SingleCellExperiment', 'scater', 'scran', 'celldex', 'SingleR', 'org.Hs.eg.db'))"

# Copy the package code into the image
COPY . /molecular
WORKDIR /molecular

# Install your package from local source
RUN R -e "devtools::install('/molecular')"

# Default command to run your analysis
CMD ["Rscript", "run_analysis.R"]
