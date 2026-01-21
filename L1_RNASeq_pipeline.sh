#!/bin/bash

################################################################################
# L1 RNA-Seq Alignment Pipeline
################################################################################
# Description: Pipeline for quantifying L1 retrotransposon and gene expression
#              from paired-end RNA-Seq data using stringent Bowtie1 alignment
#
# Author: Minh Tran
# Version: 1.0
# Last Modified: 2026-01-04
# Documentation: https://github.com/minh0620tran/L1_RNASeq_Pipeline
# Issues: https://github.com/minh0620tran/L1-RNASeq-Pipeline/issues
################################################################################

################################################################################
# SLURM Configuration - CUSTOMIZE FOR YOUR CLUSTER
################################################################################
# These settings are optimized for:
# - 128GB RAM per node (minimum 64GB)
# - 12+ CPU cores
# - Long-running jobs (12-14 hours typical, up to 48 hours)
#
# REQUIRED CHANGES before first run:
# ┌─────────────────────────────────────────────────────────────────────────┐
# │ 1. Line 46: Change YOUR_EMAIL_HERE to your email                        │
# │ 2. Line 77: Set path to your Singularity container                      │
# │ 3. Line 114: Set bind paths for your HPC filesystem                     │
# └─────────────────────────────────────────────────────────────────────────┘
#
# OPTIONAL adjustments:
# - Line 49: --time (reduce if your cluster has shorter limits)
# - Line 48: Add --partition if your HPC requires it ((e.g., 'centos7', 'gpu', 'highmem')(uncomment the line))
# - Line 51: --ntasks-per-node (adjust based on available cores)
# - Line 52: --mem (reduce if nodes have less RAM, minimum 64000)
# Common HPC variations:
# - QOS: 'long', 'standard', 'normal', 'batch', 'general'
# - Memory format: '--mem=128G' or '--mem=128000' (MB)
# - Time format: '--time=7-00:00:00' or '--time=168:00:00' (hours)
################################################################################

#SBATCH --qos=long
#SBATCH --job-name=L1-RNASeq
#SBATCH -o L1_RNASeq_OutputLog_%j.txt
#SBATCH -e L1_RNASeq_ErrorLog_%j.txt
#SBATCH --mail-user=YOUR_EMAIL_HERE #Change YOUR_EMAIL_HERE to your email
#SBATCH --mail-type=ALL
##SBATCH --partition=standard
#SBATCH --time=7-00:00:00
#SBATCH --nodes=1 
#SBATCH --ntasks-per-node=12
#SBATCH --mem=128000

################################################################################
# Pipeline Configuration
################################################################################

# Number of threads (should match --ntasks-per-node)
THREADS=12

################################################################################
# Singularity Container Setup - REQUIRED CONFIGURATION
################################################################################
# 
# STEP 1: Download the container (one-time setup)
#   singularity pull docker://minhntran/l1-rnaseq-pipeline:latest
#
# STEP 2: Set the path to your container below
# STEP 3: Configure bind mounts for your HPC filesystem
# 

# Load singularity (specify version if needed) 
module load singularity

# REQUIRED: Path to the Singularity container
# Replace with the actual path where you downloaded the .sif file
CONTAINER="/path/to/l1-rnaseq-pipeline_latest.sif"

# Verify container exists
if [ ! -f "$CONTAINER" ]; then
    echo "============================================"
    echo "ERROR: Singularity container not found"
    echo "============================================"
    echo "Expected location: $CONTAINER"
    echo ""
    echo "Please download the container first:"
    echo "  singularity pull docker://minhntran/l1-rnaseq-pipeline:latest"
    echo ""
    echo "Then update CONTAINER variable (line 82) with the correct path."
    echo "============================================"
    exit 1
fi

echo "Using Singularity container: $CONTAINER"

################################################################################
# Bind Mount Configuration
################################################################################
#
# Singularity needs to know which directories to make accessible inside 
# the container. Update the bind paths below for your HPC filesystem.
#
# Examples for common HPC setups:
#   --bind /lustre/project/yourproject:/lustre/project/yourproject
#   --bind /scratch/youruser:/scratch/youruser
#   --bind /home:/home,/data:/data
#
# Replace the path below with your actual project/data directory:
#

# Wrapper function to run tools inside container
run_tool() {
    singularity exec \
        --bind /path/to/your/project:/path/to/your/project \
        "$CONTAINER" "$@"
}

################################################################################
# Input Validation
################################################################################

# Check if correct number of arguments provided
if [ $# -ne 2 ]; then
    echo "============================================"
    echo "ERROR: Incorrect number of arguments"
    echo "============================================"
    echo "Usage: sbatch $0 <input_R1.fastq> <pipeline_directory>"
    echo ""
    echo "Example:"
    echo "  sbatch $0 sample_1.fastq.gz /path/to/L1-RNASeq-Pipeline"
    echo ""
    echo "Supported file formats:"
    echo "  - sample_1.fastq / sample_2.fastq"
    echo "  - sample_1.fastq.gz / sample_2.fastq.gz"
    echo "  - sample_R1.fastq / sample_R2.fastq"
    echo "  - sample.R1.fq.gz / sample.R2.fq.gz"
    echo "============================================"
    exit 1
fi

INPUT_FILE=$1
L1_RNASeq_dir=$2

# Check if input file exists
if [ ! -f "$INPUT_FILE" ]; then
    echo "ERROR: Input file not found: $INPUT_FILE"
    exit 1
fi

# Check if pipeline directory exists
if [ ! -d "$L1_RNASeq_dir" ]; then
    echo "ERROR: Pipeline directory not found: $L1_RNASeq_dir"
    echo "Make sure you've downloaded the reference data"
    exit 1
fi

################################################################################
# Set Reference Paths
################################################################################

# Bowtie index
hg38_bwt1_index="${L1_RNASeq_dir}/bowtie_indexes/hg38_chr_labels"

# BED files for bedtools coverage
L1_plus_bed="${L1_RNASeq_dir}/annotations/FL-L1-UCSC-plus-hg38-subfamily.bed"
L1_minus_bed="${L1_RNASeq_dir}/annotations/FL-L1-UCSC-minus-hg38-subfamily.bed"

# GTF files for featureCounts
FL_L1_PLUS_GTF="${L1_RNASeq_dir}/annotations/FL-L1-UCSC-plus-hg38-subfamily.gtf"
FL_L1_MINUS_GTF="${L1_RNASeq_dir}/annotations/FL-L1-UCSC-minus-hg38-subfamily.gtf"
HG38_GTF="${L1_RNASeq_dir}/annotations/hg38.gtf"

# Verify critical reference files exist
if [ ! -f "${hg38_bwt1_index}.1.ebwt" ]; then
    echo "ERROR: Bowtie index not found at ${hg38_bwt1_index}"
    echo "Expected files: ${hg38_bwt1_index}.*.ebwt"
    exit 1
fi

for bed_file in "$L1_plus_bed" "$L1_minus_bed"; do
    if [ ! -f "$bed_file" ]; then
        echo "ERROR: BED annotation file not found: $bed_file"
        exit 1
    fi
done

for gtf_file in "$FL_L1_PLUS_GTF" "$FL_L1_MINUS_GTF"; do
    if [ ! -f "$gtf_file" ]; then
        echo "ERROR: GTF annotation file not found: $gtf_file"
        exit 1
    fi
done

################################################################################
# Parse Input and Find Paired Files
################################################################################

# Extract sample prefix (removes _1, _2, _R1, _R2, etc.)
prefix=$(basename "$INPUT_FILE" | sed -E 's/[_.]R?[12]\.(fastq|fq)(\.gz)?$//')
fastq_file_path=$(dirname "$INPUT_FILE")

echo "=========================================="
echo "Processing sample: $prefix"
echo "=========================================="

# Find paired-end files with flexible naming
R1_file=""
R2_file=""

for suffix in "_1" "_R1" ".R1"; do
    for ext in ".fastq.gz" ".fastq" ".fq.gz" ".fq"; do
        test_r1="${fastq_file_path}/${prefix}${suffix}${ext}"
        if [ -f "$test_r1" ]; then
            R1_file="$test_r1"
            
            # Build R2 filename based on R1 suffix
            if [[ "$suffix" == "_1" ]]; then
                R2_file="${fastq_file_path}/${prefix}_2${ext}"
            elif [[ "$suffix" == "_R1" ]]; then
                R2_file="${fastq_file_path}/${prefix}_R2${ext}"
            elif [[ "$suffix" == ".R1" ]]; then
                R2_file="${fastq_file_path}/${prefix}.R2${ext}"
            fi
            
            break 2
        fi
    done
done

if [ -z "$R1_file" ] || [ ! -f "$R2_file" ]; then
    echo "ERROR: Could not find paired files for prefix: $prefix"
    echo "Searched in: $fastq_file_path"
    echo ""
    echo "Expected naming patterns:"
    echo "  ${prefix}_1.fastq / ${prefix}_2.fastq"
    echo "  ${prefix}_R1.fastq.gz / ${prefix}_R2.fastq.gz"
    echo "  ${prefix}.R1.fq / ${prefix}.R2.fq"
    exit 1
fi

echo "Found paired files:"
echo "  R1: $(basename $R1_file)"
echo "  R2: $(basename $R2_file)"

################################################################################
# Setup Output Directory
################################################################################

# Create sample-specific directory (using absolute path for Singularity)
sample_dir="$(cd "${fastq_file_path}" && pwd)/${prefix}"
mkdir -p "${sample_dir}/fastqc_reports"
mkdir -p "${sample_dir}/featureCounts"

echo "Output directory: $sample_dir"

# Move input files to sample directory
echo "Moving input files to output directory..."
mv "$R1_file" "${sample_dir}/"
mv "$R2_file" "${sample_dir}/"

# Change to sample directory for processing
cd "${sample_dir}" || exit 1

# Update file paths to local filenames
R1_file=$(basename "$R1_file")
R2_file=$(basename "$R2_file")

################################################################################
# RNASeq Alignment
################################################################################

echo ""
echo "=========================================="
echo "Starting L1 RNA-Seq Analysis Pipeline"
echo "=========================================="
echo "Sample: $prefix"
echo ""

################################################################################
# Step 1: Quality Control
################################################################################

echo ">>> [Step 1/6] Running FastQC quality control..."
run_tool fastqc "$R1_file" "$R2_file" -t $THREADS -o fastqc_reports/

if [ $? -ne 0 ]; then
    echo "ERROR: FastQC failed"
    exit 1
fi
echo "FastQC complete"

################################################################################
# Step 2: Alignment to hg38
################################################################################

echo ""
echo ">>> [Step 2/6] Aligning reads to hg38 genome with Bowtie1..."
echo "    Settings: -m 1 (unique reads only), -v 3 (3 mismatches), -X 600 (max insert), -y (tryhard)"

run_tool bash -c "bowtie -p $THREADS -m 1 -S -y -v 3 -X 600 --chunkmbs 8184 \
       ${hg38_bwt1_index} \
       -1 ${R1_file} \
       -2 ${R2_file} | \
       samtools view -@ $THREADS -hbuS - | \
       samtools sort -@ $THREADS -o ${prefix}_bowtie_hg38_sorted.bam"

if [ ! -f "${prefix}_bowtie_hg38_sorted.bam" ]; then
    echo "ERROR: Alignment failed - BAM file not created"
    exit 1
fi
echo "Alignment complete"

################################################################################
# Step 3: Remove PCR Duplicates
################################################################################

echo ""
echo ">>> [Step 3/6] Removing PCR duplicates..."

run_tool samtools rmdup ${prefix}_bowtie_hg38_sorted.bam ${prefix}_rmdup_bowtie_hg38_sorted.bam

if [ ! -f "${prefix}_rmdup_bowtie_hg38_sorted.bam" ]; then
    echo "ERROR: Duplicate removal failed"
    exit 1
fi
echo "Duplicates removed"

################################################################################
# Step 4: Strand Separation
################################################################################

echo ""
echo ">>> [Step 4/6] Separating alignments by strand..."

# Plus strand (SAM flags 83, 163)
run_tool bash -c "samtools view -@ $THREADS -h ${prefix}_rmdup_bowtie_hg38_sorted.bam | \
    awk 'substr(\$0,1,1) == \"@\" || \$2 == 83 || \$2 == 163 {print}' | \
    samtools view -@ $THREADS -bS - > ${prefix}_rmdup_bowtie_hg38_sorted_topstrand.bam"

# Minus strand (SAM flags 99, 147)
run_tool bash -c "samtools view -@ $THREADS -h ${prefix}_rmdup_bowtie_hg38_sorted.bam | \
    awk 'substr(\$0,1,1) == \"@\" || \$2 == 99 || \$2 == 147 {print}' | \
    samtools view -@ $THREADS -bS - > ${prefix}_rmdup_bowtie_hg38_sorted_bottomstrand.bam"

echo "Strand separation complete"

################################################################################
# Step 5: L1 Quantification with BEDtools
################################################################################

echo ""
echo ">>> [Step 5/6] Quantifying L1 elements with BEDtools..."

# Plus strand
run_tool bedtools coverage -a ${L1_plus_bed} \
    -b ${prefix}_rmdup_bowtie_hg38_sorted_topstrand.bam > \
    ${prefix}_rmdup_bowtie_hg38_sorted_bowtie_tryhard_plus_top.txt

echo "Plus strand: $(wc -l < ${prefix}_rmdup_bowtie_hg38_sorted_bowtie_tryhard_plus_top.txt) L1 elements"

# Minus strand
run_tool bedtools coverage -a ${L1_minus_bed} \
    -b ${prefix}_rmdup_bowtie_hg38_sorted_bottomstrand.bam > \
    ${prefix}_rmdup_bowtie_hg38_sorted_bowtie_tryhard_minus_bottom.txt

echo "Minus strand: $(wc -l < ${prefix}_rmdup_bowtie_hg38_sorted_bowtie_tryhard_minus_bottom.txt) L1 elements"

################################################################################
# Step 6: Gene/L1 Quantification with featureCounts
################################################################################

echo ""
echo ">>> [Step 6/6] Running featureCounts for gene/L1 quantification..."

# L1 plus strand
run_tool featureCounts -p --countReadPairs -B -C --primary -s 2 -t exon -g gene_id \
    -T $THREADS \
    -a "${FL_L1_PLUS_GTF}" \
    -o "featureCounts/${prefix}_FL_L1_plus_counts.txt" \
    ${prefix}_rmdup_bowtie_hg38_sorted_topstrand.bam

# L1 minus strand
run_tool featureCounts -p --countReadPairs -B -C --primary -s 2 -t exon -g gene_id \
    -T $THREADS \
    -a "${FL_L1_MINUS_GTF}" \
    -o "featureCounts/${prefix}_FL_L1_minus_counts.txt" \
    ${prefix}_rmdup_bowtie_hg38_sorted_bottomstrand.bam

# hg38 genes
run_tool featureCounts -p --countReadPairs -B -C --primary -s 2 -t exon -g gene_id \
    -T $THREADS \
    -a "${HG38_GTF}" \
    -o "featureCounts/${prefix}_hg38_genes_counts.txt" \
    ${prefix}_rmdup_bowtie_hg38_sorted.bam

echo "featureCounts complete"

################################################################################
# Index BAM Files
################################################################################

echo ""
echo ">>> Indexing BAM files..."

run_tool samtools index ${prefix}_rmdup_bowtie_hg38_sorted.bam

echo "BAM files indexed"

################################################################################
# Completion
################################################################################

echo ""
echo "=========================================="
echo "Pipeline completed successfully!"
echo "=========================================="
echo "Sample: $prefix"
echo "Output directory: ${sample_dir}"
echo "Runtime: $SECONDS seconds"
echo ""
echo "Output files:"
echo ""
echo "Alignments:"
echo "  • ${prefix}_rmdup_bowtie_hg38_sorted.bam"
echo ""
echo "L1 Quantification (for L1ABA):"
echo "  • ${prefix}_rmdup_bowtie_hg38_sorted_bowtie_tryhard_plus_top.txt"
echo "  • ${prefix}_rmdup_bowtie_hg38_sorted_bowtie_tryhard_minus_bottom.txt"
echo "  • ${prefix}_rmdup_bowtie_hg38_sorted_topstrand.bam"
echo "  • ${prefix}_rmdup_bowtie_hg38_sorted_bottomstrand.bam"
echo ""
echo "Gene/L1 Counts (featureCounts):"
echo "  • featureCounts/${prefix}_FL_L1_plus_counts.txt"
echo "  • featureCounts/${prefix}_FL_L1_minus_counts.txt"
echo "  • featureCounts/${prefix}_hg38_genes_counts.txt"
echo ""
echo "Quality Control:"
echo "  • fastqc_reports/"
echo ""
echo "=========================================="