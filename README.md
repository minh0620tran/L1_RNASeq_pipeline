# L1 RNA-Seq Alignment Pipeline


Pipeline for quantifying LINE-1 (L1) retrotransposon and gene expression from stranded, paired-end RNA-Seq data using stringent Bowtie1 alignment.

## Overview

Due to the highly repetitive nature of LINE-1 elements, this pipeline uses **Bowtie1** with unique alignment settings (`-m 1`) to ensure high specificity in mapping reads to L1 loci. The pipeline is optimized for High-Performance Computing (HPC) environments and uses **Singularity containers** to ensure reproducibility across different systems.

The outputs of this alignment will then be used for downstream L1 Automated Bioinformatics Analysis (L1ABA)

### Pipeline Steps

1. **Quality Control** - FastQC for read quality assessment
2. **Alignment** - Bowtie1 with stringent unique mapping to hg38
3. **Duplicate Removal** - Remove PCR duplicates with samtools
4. **Strand Separation** - Separate plus/minus strand alignments
5. **L1 Quantification** - BEDtools coverage for L1 elements (for L1ABA downstream analysis)
6. **Gene Quantification** - featureCounts for genes and L1 elements

### Runtime & Resources

- **Typical runtime:** 12-14 hours per sample (up to 48 hours for large samples)
- **Memory required:** 128GB RAM (minimum 64GB)
- **CPU cores:** 12 threads recommended
- **Disk space:** ~50GB per sample

---

## Installation

### Prerequisites

- HPC cluster with SLURM scheduler
- Singularity (v3.0+) - usually available as a module
- Minimum 64GB RAM per node (128GB recommended)
- ~500GB disk space for reference data and outputs

### Step 1: Clone Repository

```bash
git clone https://github.com/minh0620tran/L1_RNASeq_Pipeline.git
cd L1_RNASeq_Pipeline
```

### Step 2: Download Singularity Container

```bash
# Download pre-built container from Docker Hub
module load singularity
singularity pull docker://minhntran/l1-rnaseq-pipeline:latest

# This creates: l1-rnaseq-pipeline_latest.sif
```

The container includes all required tools:
- FastQC
- Bowtie v1.3.1
- Samtools
- BEDtools
- Subread (featureCounts)

### Step 3: Configure Bind Mounts (IMPORTANT!)

Singularity needs to know which directories to make accessible inside the container. 

**Determine your project directory:**
````bash
# Where is your data located?
pwd
# Example output: /lustre/project/yourlab/rnaseq

# Your bind path should include this directory
# For example: /lustre/project/yourlab
````

**Test that bind mounts work:**
````bash
# Test if container can see your files
singularity exec --bind /lustre/project/yourlab:/lustre/project/yourlab \
  l1-rnaseq-pipeline_latest.sif \
  ls $(pwd)
````

**Common bind mount patterns:**
- `/lustre/project/yourlab:/lustre/project/yourlab` 
- `/scratch/username:/scratch/username` (if using scratch)
- `/home:/home,/data:/data` (multiple directories, comma-separated)

**Save your bind mount command for later steps (Step 4-6)!**

### Step 4: Build hg38 Bowtie Index

Now that bind mounts are configured, build the index:
````bash
cd bowtie_indexes/

# Download hg38 reference genome
wget http://hgdownload.cse.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz
gunzip hg38.fa.gz

# Build Bowtie1 index using your configured bind mount
# Replace /lustre/project/yourlab with YOUR bind path from Step 3
singularity exec --bind /lustre/project/yourlab:/lustre/project/yourlab \
  ../l1-rnaseq-pipeline_latest.sif \
  bowtie-build hg38.fa hg38_chr_labels
````

### Step 5: Download hg38 GTF file

```bash
cd ../annotations/

wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_49/gencode.v49.annotation.gtf.gz
gunzip gencode.v49.annotation.gtf.gz
mv gencode.v49.annotation.gtf hg38.gtf
```

**Note:** L1 annotation files (BED/GTF) are included in the repository under /annotation.

### Step 6: Configure the Pipeline Script

Edit `L1_RNASeq_pipeline.sh` and update these lines:

```bash
# Line 46: Your email address
#SBATCH --mail-user=YOUR_EMAIL_HERE

# Line 48: Add partition if your HPC requires it (uncomment and modify)
# #SBATCH --partition=standard

# Line 77: Path to Singularity container
CONTAINER="/path/to/l1-rnaseq-pipeline_latest.sif"

# Line 114: Bind mount for your HPC filesystem
# Replace with your actual project/data directory
--bind /path/to/your/project:/path/to/your/project
```

**Example bind paths for common HPC setups:**
- `/lustre/project/yourlab:/lustre/project/yourlab`
- `/scratch/yourusername:/scratch/yourusername`
- `/home:/home` (bind entire home directory)

---

## Usage

### Input Requirements

**File format:**
- Paired-end FASTQ files (gzipped or uncompressed)
- Supported extensions: `.fastq`, `.fastq.gz`, `.fq`, `.fq.gz`

**Naming conventions (any of these work):**
- `sample_1.fastq` / `sample_2.fastq`
- `sample_1.fastq.gz` / `sample_2.fastq.gz`
- `sample_R1.fastq` / `sample_R2.fastq`
- `sample.R1.fq.gz` / `sample.R2.fq.gz`

### Running a Single Sample

```bash
sbatch L1_RNASeq_pipeline.sh sample_1.fastq.gz /path/to/L1_RNASeq_Pipeline
```

Arguments:
1. Path to R1 FASTQ file (R2 is auto-detected)
2. Path to pipeline directory (where reference data is located)

### Running Multiple Samples

```bash
# Process all samples with _1.fastq.gz suffix
for r1 in *_1.fastq.gz; do
  sbatch L1_RNASeq_pipeline.sh $r1 /path/to/L1_RNASeq_Pipeline
done

# Or with custom job names
for r1 in *_1.fastq.gz; do
  sample=$(basename $r1 _1.fastq.gz)
  sbatch --job-name=$sample L1_RNASeq_pipeline.sh $r1 /path/to/L1_RNASeq_Pipeline
done
```

## Output Files

For each sample, the pipeline creates a directory with the following structure:

```
sample_name/
├── fastqc_reports/                   # Quality control reports
│   ├── sample_1_fastqc.html
│   └── sample_2_fastqc.html
├── featureCounts/                    # Gene/L1 quantification
│   ├── sample_FL_L1_plus_counts.txt
│   ├── sample_FL_L1_minus_counts.txt
│   └── sample_hg38_genes_counts.txt
├── sample_rmdup_bowtie_hg38_sorted.bam           
├── sample_rmdup_bowtie_hg38_sorted.bam.bai       
├── sample_rmdup_bowtie_hg38_sorted_topstrand.bam                         #Strand-separated alignment file for topstrand
├── sample_rmdup_bowtie_hg38_sorted_bottomstrand.bam                      #Strand-separated alignment file for bottomstrand
├── sample_rmdup_bowtie_hg38_sorted_bowtie_tryhard_plus_top.txt           # L1 counts (for downstream L1ABA)
├── sample_rmdup_bowtie_hg38_sorted_bowtie_tryhard_minus_bottom.txt       # L1 counts (for downstream L1ABA)
├── sample_1.fastq.gz                             # Input files 
└── sample_2.fastq.gz
```

---

## Pipeline Details

### Alignment Strategy

The pipeline uses Bowtie1 with the following parameters optimized for L1 elements:

```bash
bowtie -p 12 -m 1 -S -y -v 3 -X 600
```

- `-m 1` : Report only **uniquely mapping reads** (critical for repetitive elements)
- `-v 3` : Allow up to 3 mismatches
- `-X 600` : Maximum insert size 600bp
- `-y` : Exhaustive search for valid alignments

### Strand-Specific Processing

1. Separates alignments by strand (SAM flags 83/163 vs 99/147)
2. Quantifies L1 elements separately for each strand
3. Matches strand-specific L1 annotations

---

## Customization

### Adjusting SLURM Parameters

If your HPC has different resource limits, adjust these lines in the script:

```bash
#SBATCH --qos=long              # Change to your QOS name
#SBATCH --time=7-00:00:00       # Reduce if needed (minimum 1 day recommended)
#SBATCH --mem=128000            # Reduce to 64000 if needed (minimum)
#SBATCH --ntasks-per-node=12    # Adjust based on available cores
```

### Memory Requirements by Sample Size

| FASTQ Size (Gzipped) | Recommended RAM | Typical Runtime |
|---------------------|-----------------|-----------------|
| < 5 GB per pair     | 64 GB           | 8-10 hours      |
| 5-10 GB per pair    | 128 GB          | 12-14 hours     |
| > 10 GB per pair    | 256 GB          | 24-48 hours     |

---

### Getting Help

- **Issues:** https://github.com/minh0620tran/L1-RNASeq-Pipeline/issues
- **Email:**  mtran13@tulane.edu

---

## Citation

If you use this pipeline in your research, please cite:

```
Deininger P, Morales ME, White TB, Baddoo M, Hedges DJ, Servant G, Srivastav S, Smither ME, Concha M, DeHaro DL, Flemington EK, Belancio VP. A comprehensive approach to expression of L1 loci. Nucleic Acids Res. 2017 Mar 17;45(5):e31. doi: 10.1093/nar/gkw1067. PMID: 27899577; PMCID: PMC5389711.
Kaul T, Morales ME, Smither E, Baddoo M, Belancio VP, Deininger P. RNA Next-Generation Sequencing and a Bioinformatics Pipeline to Identify Expressed LINE-1s at the Locus-Specific Level. J Vis Exp. 2019 May 19;(147):10.3791/59771. doi: 10.3791/59771. PMID: 31157783; PMCID: PMC7371004.

```

---

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.
