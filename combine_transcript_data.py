import csv
import argparse
import re
from Bio import SeqIO
from Bio.Seq import Seq

def find_longest_orf(sequence):
    longest_orf = ''
    for frame in range(3):
        translated = Seq(sequence[frame:]).translate(cds=False)
        orfs = re.findall(r'(?=(M[^*]*))', str(translated))
        for orf in orfs:
            if len(orf) > len(longest_orf):
                longest_orf = orf

    for frame in range(3):  # Iterate through three reading frames for the reverse strand
        rev_translated = Seq(sequence[:len(sequence)-frame]).reverse_complement().translate(cds=False)
        rev_orfs = re.findall(r'(?=(M[^*]*))', str(rev_translated))
        for orf in rev_orfs:
            if len(orf) > len(longest_orf):
                longest_orf = orf

    return longest_orf

def combine_data(fasta_file, quant_file, output_file):
    # Read the fasta file
    sequences = {}
    for record in SeqIO.parse(fasta_file, "fasta"):
        sequences[record.id] = str(record.seq)

    # Read the quant.sf file and create the output CSV file
    with open(quant_file) as quant, open(output_file, 'w') as output:
        csv_reader = csv.DictReader(quant, delimiter='\t')
        csv_writer = csv.writer(output)
        csv_writer.writerow(['transcript_id', 'sequence', 'length', 'effectivelength', 'tpm', 'numreads', 'longest_coding_orf'])

        for row in csv_reader:
            transcript_id = row['Name']
            if transcript_id in sequences:
                sequence = sequences[transcript_id]
                longest_orf = find_longest_orf(sequence)
                csv_writer.writerow([transcript_id, sequence, row['Length'], row['EffectiveLength'], row['TPM'], row['NumReads'], longest_orf])
            
def main():
    parser = argparse.ArgumentParser(description='Combine fasta sequences and quant.sf data into a CSV file.')
    parser.add_argument('fasta_file', help='Input fasta file containing transcript sequences.')
    parser.add_argument('quant_file', help='Input quant.sf file containing expression data.')
    parser.add_argument('output_file', help='Output CSV file to store the combined data.')

    args = parser.parse_args()

    combine_data(args.fasta_file, args.quant_file, args.output_file)

if __name__ == '__main__':
    main()
