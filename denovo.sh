# work in progress - the STAR command needs way too much ram
# will try --genomeSAsparseD plus filtering beforehand

mkdir -p ./data/spades_output
mkdir -p ./STAR_index/rnaspades_nogenomeref
mkdir -p ./data/aligned_reads/rnaspades_nogenomeref

rnaspades.py \
--pe1-1 data/fastq/SRR2273875_1.fastq --pe1-2 data/fastq/SRR2273875_2.fastq \
-o ./data/spades_output

STAR \
--runMode genomeGenerate \
--limitGenomeGenerateRAM 51000000000 \
--genomeDir ./STAR_index/rnaspades_nogenomeref \
--genomeFastaFiles ./data/spades_output/transcripts.fasta \
--runThreadN 12


STAR \
--genomeDir ./STAR_index/rnaspades_nogenomeref \
--readFilesIn ./data/fastq/SRR2273875_1.fastq ./data/fastq/SRR2273875_2.fastq \
--runThreadN numberOfThreads 12 \
--outFileNamePrefix ./data/aligned_reads/rnaspades_nogenomeref