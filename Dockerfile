FROM rocker/r-base

WORKDIR /home/molecular

COPY final.R .

RUN R -e "install.packages(c('ggplot2', 'stringr', 'igraph', 'BiocManager'), repos='https://cloud.r-project.org')"

RUN R -e "BiocManager::install(c('SingleCellExperiment', 'DropletUtils', 'scater', 'scran', 'org.Hs.eg.db', 'celldex', 'SingleR'))"

CMD ["Rscript", "final.R"]
