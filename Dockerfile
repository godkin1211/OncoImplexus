# 1. Base Image - R 4.4 with Shiny and Tidyverse
FROM rocker/shiny-verse:4.4.0

# 2. System Dependencies
RUN apt-get update && apt-get install -y \
    libxml2-dev \
    libcurl4-openssl-dev \
    libssl-dev \
    zlib1g-dev \
    libbz2-dev \
    liblzma-dev \
    libglpk-dev \
    libpq-dev \
    libgsl-dev \
    libcairo2-dev \
    libxt-dev \
    && rm -rf /var/lib/apt/lists/*

# 3. Bioconductor Setup
RUN R -e "install.packages('BiocManager', repos='http://cran.rstudio.com/')"

# 4. Install Heavy Genome Dependencies (Cached Layer)
# Installing hg19, hg38, and annotation packages. This step takes a long time.
RUN apt-get update && apt-get install -y \
    libssl-dev \
    libxml2-dev \
    libcurl4-openssl-dev \
    libgit2-dev \
    libglpk-dev \
    libcairo2-dev \
    libxt-dev \
    libpq-dev \
    libgsl-dev \
    libbz2-dev \
    liblzma-dev \
    libfftw3-dev \
    libopenblas-dev \
    libtiff5-dev \
    libjpeg-dev \
    libpng-dev \
    libssh2-1-dev \
    zlib1g-dev \
    && rm -rf /var/lib/apt/lists/*

RUN R -e "BiocManager::install(c( \
    'GenomicRanges', \
    'StructuralVariantAnnotation', \
    'BSgenome', \
    'BSgenome.Hsapiens.UCSC.hg19', \
    'BSgenome.Hsapiens.UCSC.hg38', \
    'TxDb.Hsapiens.UCSC.hg19.knownGene', \
    'TxDb.Hsapiens.UCSC.hg38.knownGene', \
    'org.Hs.eg.db', \
    'VariantAnnotation', \
    'graph', \
    'RBGL' \
    ))"

# 5. Install other CRAN dependencies for OncoImplexus
RUN R -e "install.packages(c('shinydashboard', 'DT', 'colourpicker', 'gridExtra', 'circlize', 'cowplot'))"

# 6. Install OncoImplexus Package
# Copy the entire project into the container
COPY . /app/OncoImplexus

# Install the package from source
RUN R -e "install.packages('matrixStats', repos='http://cran.rstudio.com/')" \
    && R -e "BiocManager::install(c('GenomicFeatures', 'StructuralVariantAnnotation', 'BSgenome', 'VariantAnnotation'))" \
    && R -e "devtools::install('/app/OncoImplexus', dependencies=TRUE)"

# 7. Setup Shiny App
# Clean default app
RUN rm -rf /srv/shiny-server/*

# Copy the Shiny app file
COPY inst/shiny/app.R /srv/shiny-server/

# Copy Gene Annotation RDS files explicitly to the app directory
# This ensures the app can find them in its working directory
COPY hg19_genes.rds /srv/shiny-server/
COPY hg38_genes.rds /srv/shiny-server/

# 8. Exposure and Command
EXPOSE 3838
CMD ["/usr/bin/shiny-server"]
