# Bowtie Index Files

This directory will contain the hg38 Bowtie1 index files after you build them.

## Build Instructions

See the main [README](../README.md#step-4-build-hg38-bowtie-index) for complete setup instructions.

**Quick summary:**
```bash
cd bowtie_indexes/

# Download reference
wget http://hgdownload.cse.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz
gunzip hg38.fa.gz

# Build index (configure bind mounts first - see main README!)
singularity exec --bind /your/project:/your/project \
  ../l1-rnaseq-pipeline_latest.sif \
  bowtie-build hg38.fa hg38_chr_labels
```

## Expected Files

After building, you should have 6 files (~4GB total):
- `hg38_chr_labels.{1-4}.ebwt`
- `hg38_chr_labels.rev.{1,2}.ebwt`

**Note:** These files are not included in the repository due to size. Each user must build them once.