#!/bin/bash

# stop if there is an error
set -e

# verbose output
set -x
set -v

# dependency check
for lib in STAR samtools prefetch fasterq-dump stringtie; do
  if ! command -v "$lib" &> /dev/null; then
    echo "Error: $lib is required but not found. Please install $lib and ensure it's in your PATH."
    exit 1
  fi
done

# Set up directories
DATA_DIR="./data"
SRA_DIR="$DATA_DIR/sra_dump"
FASTQ_DIR="$DATA_DIR/fastq"
STAR_INDEX="./STAR_index"
ALIGNED_DIR="$DATA_DIR/aligned_reads"

for dir in {$DATA_DIR,$SRA_DIR,$FASTQ_DIR,$STAR_INDEX,$ALIGNED_DIR}; do
  mkdir -p $dir
done

# Get reference genome and annotation
GENOME_FASTA="GCF_001888835.1_ASM188883v1_genomic.fna"
GENOME_GTF="GCF_001888835.1_ASM188883v1_genomic.gtf"

if [ ! -f "$GENOME_FASTA" ]; then
    wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/annotation_releases/89399/100/GCF_001888835.1_ASM188883v1/GCF_001888835.1_ASM188883v1_genomic.fna.gz
    gunzip GCF_001888835.1_ASM188883v1_genomic.fna.gz
fi

if [ ! -f "$GENOME_GTF" ]; then
    wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/annotation_releases/89399/100/GCF_001888835.1_ASM188883v1/GCF_001888835.1_ASM188883v1_genomic.gtf.gz
    gunzip GCF_001888835.1_ASM188883v1_genomic.gtf.gz
fi

# Index genome - watch the RAM allocation (Needs at least 50 gb)
if [ ! -f "$STAR_INDEX/SAindex" ]; then
    STAR --runMode genomeGenerate \
         --limitGenomeGenerateRAM 51000000000 \
         --genomeDir $STAR_INDEX \
         --genomeFastaFiles $GENOME_FASTA \
         --runThreadN 8 \
         --sjdbGTFfile $GENOME_GTF \
         --sjdbOverhang 100
fi

# RNA seq retrieval and processing
# (note the last two are from R. ferrumequinum, but these were actually used for
# what you see on ncbi)
RUNS=(
  "SRR1048142"
  "SRR1048140"
  "SRR1584445"
  "SRR1584446"
  "SRR1584447"
  "SRR2153217"
  "SRR2273738"
  "SRR2273739"
  "SRR2273740"
  "SRR2273762"
  "SRR2273816"
  "SRR2273875"
  "SRR2273931"
  "SRR2754983"
  "SRR2757329"
)

# download
for srr in "${RUNS[@]}"
do
  if [ ! -f "$SRA_DIR/$srr.sra" ]; then
    prefetch -O "$SRA_DIR" "$srr"
  fi
done

# convert to fastq
for sra_file in $SRA_DIR/*/
do
  sra_base=$(basename "${sra_file%/}")
  if [ ! -f "$FASTQ_DIR/${sra_base}_1.fastq" ] || [ ! -f "$FASTQ_DIR/${sra_base}_2.fastq" ]; then
    fasterq-dump -v -O $FASTQ_DIR --split-3 "$sra_file"
  fi
done

######## Missing some QC! Should run at least a MultiQC report to check ########

# STAR alignments against reference genome - outputs aligned reads in bam format
for read1 in $FASTQ_DIR/*_1.fastq; do
    read2="${read1/_1.fastq/_2.fastq}"
    prefix=$(basename "${read1%_1.fastq}")
    
    if [ ! -f "$ALIGNED_DIR/${prefix}_Aligned.sortedByCoord.out.bam" ]; then
        STAR --genomeDir $STAR_INDEX \
             --readFilesIn $read1 $read2 \
             --runThreadN 8 \
             --outSAMtype BAM SortedByCoordinate \
             --outFileNamePrefix $ALIGNED_DIR/${prefix}_ \
             --outFilterType BySJout \
             --alignSJoverhangMin 8
    fi
done

# you can preview them with samtools view -h filename.bam | head -n 10

# Merge BAMs. Probably need to be more thoughtful about just merging all files.
# It might make more sense to group runs from the same lab/sequencer.
# (plus for ifitm, SRR2273740 and SRR2273875 may be sufficient.)
if [ ! -f "$ALIGNED_DIR/merged.bam" ]; then
    samtools merge -f $ALIGNED_DIR/merged.bam $ALIGNED_DIR/*_Aligned.sortedByCoord.out.bam
fi

# Index BAMs
for bam_file in $ALIGNED_DIR/*.bam; do
    if [ ! -f "${bam_file}.bai" ]; then
        samtools index $bam_file
    fi
done

# Assemble transcripts - uses reference genome annotation file
for bam_file in $ALIGNED_DIR/*.bam; do
    sample_name=$(basename "${bam_file%.*}")
    if [ ! -f "$ALIGNED_DIR/${sample_name}.gtf" ]; then
        stringtie -p 8 -G $GENOME_GTF -o $ALIGNED_DIR/${sample_name}.gtf -l ${sample_name} $bam_file
    fi
done