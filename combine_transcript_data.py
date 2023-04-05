import csv
import argparse
from Bio import SeqIO

def combine_data(fasta_file, quant_file, output_file):
    # Read the fasta file
    sequences = {}
    for record in SeqIO.parse(fasta_file, "fasta"):
        sequences[record.id] = str(record.seq)

    # Read the quant.sf file and create the output CSV file
    with open(quant_file) as quant, open(output_file, 'w') as output:
        csv_reader = csv.DictReader(quant, delimiter='\t')
        csv_writer = csv.writer(output)
        csv_writer.writerow(['transcript_id', 'sequence', 'length', 'effectivelength', 'tpm', 'numreads'])

        for row in csv_reader:
            transcript_id = row['Name']
            if transcript_id in sequences:
                csv_writer.writerow([transcript_id, sequences[transcript_id], row['Length'], row['EffectiveLength'], row['TPM'], row['NumReads']])
            
def main():
    parser = argparse.ArgumentParser(description='Combine fasta sequences and quant.sf data into a CSV file.')
    parser.add_argument('fasta_file', help='Input fasta file containing transcript sequences.')
    parser.add_argument('quant_file', help='Input quant.sf file containing expression data.')
    parser.add_argument('output_file', help='Output CSV file to store the combined data.')

    args = parser.parse_args()

    combine_data(args.fasta_file, args.quant_file, args.output_file)

if __name__ == '__main__':
    main()
