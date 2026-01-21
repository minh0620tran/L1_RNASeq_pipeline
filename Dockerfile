# L1 RNA-seq Pipeline Container
FROM continuumio/miniconda3:24.1.2-0

# Install system dependencies
RUN apt-get update && apt-get install -y \
    wget \
    unzip \
    && rm -rf /var/lib/apt/lists/*

# Copy environment file
COPY environment.yml /tmp/environment.yml

# Create conda environment
RUN conda env create -f /tmp/environment.yml && \
    conda clean --all --yes

# Put conda env in PATH
ENV PATH=/opt/conda/envs/l1-rnaseq-pipeline/bin:$PATH

# Set working directory
WORKDIR /data

# Test that tools work 
RUN fastqc --version && \
    bowtie --version && \
    samtools --version && \
    bedtools --version && \
    featureCounts -v

CMD ["/bin/bash"]