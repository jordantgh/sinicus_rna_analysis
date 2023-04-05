#!/bin/bash

# stop if there is an error
set -e

# verbose output
set -x
set -v

SALMON_INDEX="./salmon_index/rnaspades_nogenomeref"
DATA_DIR="./data"
ALIGNED_DIR="$DATA_DIR/aligned_reads"
SPADES_OUTPUT="$DATA_DIR/spades_output"
SALMON_OUTPUT="$DATA_DIR/aligned_reads/salmon_output"
FINAL_OUTPUT="$DATA_DIR/denovo_results"
FASTQ_DIR="$DATA_DIR/fastq"

for dir in {$SPADES_OUTPUT,$SALMON_INDEX,$ALIGNED_DIR,$SALMON_OUTPUT,$FINAL_OUTPUT}; do
  mkdir -p $dir
done

# Define the runs for processing
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

# RNA seq processing
for run in "${RUNS[@]}"; do
  # Run rnaSPAdes
  if [ ! -d "$SPADES_OUTPUT/$run" ]; then
    rnaspades.py \
      --pe1-1 $FASTQ_DIR/${run}_1.fastq \
      --pe1-2 $FASTQ_DIR/${run}_2.fastq \
      -o $SPADES_OUTPUT/$run
  fi

  # Index with Salmon
  if [ ! -d "$SALMON_INDEX/$run" ]; then
    salmon index -t $SPADES_OUTPUT/$run/hard_filtered_transcripts.fasta \
    -i $SALMON_INDEX/$run
  fi

  # Quantify with Salmon
  if [ ! -d "$SALMON_OUTPUT/$run" ]; then
    salmon quant -i $SALMON_INDEX/$run -l A \
      -1 $FASTQ_DIR/${run}_1.fastq -2 $FASTQ_DIR/${run}_2.fastq \
      -p 12 -o $SALMON_OUTPUT/$run
  fi

  # Create BLAST database of transcripts
  if [ ! -d "$FINAL_OUTPUT/$run/blast_db" ]; then
    mkdir -p $FINAL_OUTPUT/$run/blast_db
    makeblastdb -in $SPADES_OUTPUT/$run/hard_filtered_transcripts.fasta \
    -dbtype nucl -out $FINAL_OUTPUT/$run/blast_db/${run}_transcripts_db
  fi

  # Perform blast search
  QUERY_FILE="sinicus_pred_ifitm3_x2.fasta"
  if [ ! -f "$QUERY_FILE" ]; then
    echo "Error: Query file not found: $QUERY_FILE"
    exit 1
  fi
  
  if [ ! -f "$FINAL_OUTPUT/$run/${run}_blast_results.txt" ]; then
    blastn -query $QUERY_FILE \
           -db $FINAL_OUTPUT/$run/blast_db/${run}_transcripts_db \
           -out $FINAL_OUTPUT/$run/${run}_blast_results.txt \
           -outfmt 6 \
           -evalue 1e-5 \
           -max_target_seqs 10
  fi

  # Extract transcript IDs from blast results
  if [ ! -f "$FINAL_OUTPUT/$run/${run}_transcript_ids.txt" ]; then
    cut -f 2 $FINAL_OUTPUT/$run/${run}_blast_results.txt \
    > $FINAL_OUTPUT/$run/${run}_transcript_ids.txt
  fi
  
  # Filter the original fasta file with the extracted transcript IDs
  if [ ! -f "$FINAL_OUTPUT/$run/${run}_subset_transcripts.fasta" ]; then
    seqtk subseq \
    $SPADES_OUTPUT/$run/hard_filtered_transcripts.fasta \
    $FINAL_OUTPUT/$run/${run}_transcript_ids.txt \
    > $FINAL_OUTPUT/$run/${run}_subset_transcripts.fasta
  fi
  
  # Create a table with transcript id, sequence, and quant.sf data (length,   effectivelength, tpm, numreads)
  if [ ! -f "$FINAL_OUTPUT/$run/${run}_combined_data.csv" ]; then
    poetry run python ./aux_scripts/tabulate_tx.py \
    $FINAL_OUTPUT/$run/${run}_subset_transcripts.fasta \
    $SALMON_OUTPUT/$run/quant.sf $FINAL_OUTPUT/$run/${run}_combined_data.csv
  fi
done