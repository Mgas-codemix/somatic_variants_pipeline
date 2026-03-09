# Somatic Variants Pipeline - reproducible image with all tools
# Build: docker build -t somatic-variants:latest -f docker/Dockerfile .
# Optional: place MuTect JAR in docker/jars/ and ANNOVAR tarball as docker/annovar.tar.gz for full from_fastq support

FROM ubuntu:22.04

LABEL org.opencontainers.image.title="somatic-variants-pipeline" \
      org.opencontainers.image.description="Somatic Variants Pipeline — bioinformatics tools and R packages" \
      org.opencontainers.image.url="https://github.com/Mgas-codemix/somatic_variants_pipeline" \
      org.opencontainers.image.licenses="MIT"

ENV DEBIAN_FRONTEND=noninteractive
ENV LANG=C.UTF-8

# Java 8 for GATK 3, Picard, MuTect, VarScan
RUN apt-get update && apt-get install -y --no-install-recommends \
    openjdk-8-jre-headless \
    ca-certificates \
    wget \
    curl \
    gzip \
    make \
    gcc \
    libz-dev \
    libncurses5-dev \
    libncurses-dev \
    libbz2-dev \
    liblzma-dev \
    libcurl4-openssl-dev \
    libssl-dev \
    perl \
    python3 \
    python3-pip \
    r-base \
    r-base-dev \
    && rm -rf /var/lib/apt/lists/*

# BWA, Samtools, Cutadapt, FastQC from apt (Ubuntu 22.04)
RUN apt-get update && apt-get install -y --no-install-recommends \
    bwa \
    samtools \
    cutadapt \
    fastqc \
    && rm -rf /var/lib/apt/lists/*

# seqtk (not always in apt)
RUN wget -q https://github.com/lh3/seqtk/archive/refs/tags/v1.4.tar.gz -O /tmp/seqtk.tar.gz \
    && tar -xzf /tmp/seqtk.tar.gz -C /tmp \
    && make -C /tmp/seqtk-1.4 \
    && mv /tmp/seqtk-1.4/seqtk /usr/local/bin/ \
    && rm -rf /tmp/seqtk*

# R packages required by pipeline scripts
RUN R -e "install.packages(c('data.table', 'ggplot2', 'ggsci', 'ggpubr', 'reshape2', 'openxlsx', 'pROC'), repos='https://cloud.r-project.org/', Ncpus=2)"

# Bioconductor packages for real BED interval overlap
RUN R -e "if (!requireNamespace('BiocManager', quietly=TRUE)) install.packages('BiocManager', repos='https://cloud.r-project.org/'); BiocManager::install(c('GenomicRanges', 'IRanges', 'S4Vectors'), ask=FALSE)"

# MultiQC (Python-based QC aggregator)
RUN pip3 install multiqc

# BCFtools
RUN apt-get update && apt-get install -y --no-install-recommends bcftools && rm -rf /var/lib/apt/lists/*

# Quarto (for HTML clinical report rendering)
RUN wget -q https://github.com/quarto-dev/quarto-cli/releases/download/v1.4.553/quarto-1.4.553-linux-amd64.deb -O /tmp/quarto.deb \
    && dpkg -i /tmp/quarto.deb && rm /tmp/quarto.deb

# R packages for Quarto clinical report
RUN R -e "install.packages(c('dplyr','tidyverse','readr','readxl','gt','htmltools','tibble','gtExtras','knitr','kableExtra','DT','tidyr','rmarkdown','ggplot2'), repos='https://cloud.r-project.org/', Ncpus=2)"

# Install Java tools to /opt
ENV OPT=/opt/somatic
RUN mkdir -p ${OPT}/jars

# Picard (Apache 2.0)
RUN wget -q https://github.com/broadinstitute/picard/releases/download/2.27.5/picard.jar -O ${OPT}/jars/picard.jar

# GATK 3.8-1 (Broad license - research use)
RUN wget -q https://storage.googleapis.com/gatk-software/package-archive/gatk/GenomeAnalysisTK-3.8-1-0-gf15c1c3ef.tar.bz2 -O /tmp/gatk.tar.bz2 \
    && tar -xjf /tmp/gatk.tar.bz2 -C /tmp \
    && mv /tmp/GenomeAnalysisTK-3.8-1-0-gf15c1c3ef/GenomeAnalysisTK.jar ${OPT}/jars/ \
    && rm -rf /tmp/gatk.tar.bz2 /tmp/GenomeAnalysisTK-*

# VarScan 2.4.6 (GitHub release; v2.4.4 assets no longer available)
RUN wget -q https://github.com/dkoboldt/varscan/releases/download/v2.4.6/VarScan.v2.4.6.jar -O ${OPT}/jars/VarScan.v2.4.6.jar

# MuTect 1.1.7: add jar to docker/jars/ (e.g. muTect-1.1.7.jar) and rebuild, or mount at run time.
# Broad license; obtain from Broad Institute or build from source.
COPY docker/jars/ ${OPT}/jars/

# ANNOVAR: add annovar.tar.gz to docker/ and uncomment below, or mount /opt/annovar at run time.
# Download from https://www.openbioinformatics.org/annovar/ (registration required).
# RUN mkdir -p /opt/annovar && tar -xzf /tmp/annovar.tar.gz -C /opt/annovar --strip-components=1
RUN mkdir -p /opt/annovar/humandb

# Use Java 8 for GATK 3 / Picard / MuTect (path is java-8-openjdk-amd64 or java-8-openjdk-arm64)
RUN /bin/bash -c 'J8=$(ls -d /usr/lib/jvm/java-8-openjdk-* 2>/dev/null | head -1) && ln -sf "$J8" /usr/lib/jvm/java-8'
ENV JAVA_HOME=/usr/lib/jvm/java-8
ENV PATH="${JAVA_HOME}/bin:${OPT}/jars:/opt/annovar:${PATH}"

# Default tool paths (Nextflow params can override)
ENV GATK_JAR=${OPT}/jars/GenomeAnalysisTK.jar
ENV PICARD_JAR=${OPT}/jars/picard.jar
ENV MUTECT_JAR=${OPT}/jars/muTect-1.1.7.jar
ENV VARSCAN_JAR=${OPT}/jars/VarScan.v2.4.6.jar
ENV ANNOVAR_BASEDIR=/opt/annovar

WORKDIR /work
