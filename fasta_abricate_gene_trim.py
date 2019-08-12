#!/bin/env python

import argparse
import csv
import sys

import Bio.SeqIO

#from pprint import pprint


def check_positive(value):
    ivalue = int(value)
    if ivalue < 0:
        raise argparse.ArgumentTypeError("%s is an invalid positive int value" % value)
    return ivalue

def main():    
    parser = argparse.ArgumentParser()
    parser.add_argument('--fasta', dest='fasta', required=True, help='input fasta')
    parser.add_argument( '--abricate-output', dest='abricate_output', required=True, help='Output file from abricate (.tsv)' )
    parser.add_argument('--buffer-size', dest='buffer_size', default=0, type=check_positive, help='Number of bases to include on either side of target gene')
    parser.add_argument('--target-gene', dest='target_gene', required=True, help='Gene to extract from fasta')
    args = parser.parse_args()

    target_gene = {}
    
    with open(args.abricate_output) as csvfile:
        abricate_output_fieldnames = [
            'file',
            'sequence',
            'start',
            'end',
            'gene',
            'coverage',
            'coverage_map',
            'gaps',
            'percent_coverage',
            'percent_identity',
            'database',
            'accession',
            'product',
        ]
        
        reader = csv.DictReader(csvfile, fieldnames=abricate_output_fieldnames, delimiter='\t')
        for row in reader:
            if row['gene'] == args.target_gene:
                target_gene = row
                target_gene['start'] = int(target_gene['start'])
                target_gene['end'] = int(target_gene['end'])
                target_gene['percent_coverage'] = float(target_gene['percent_coverage'])
                target_gene['percent_identity'] = float(target_gene['percent_identity'])

        
    for record in Bio.SeqIO.parse(args.fasta, "fasta"):
        if str(record.id) == target_gene['sequence']:
            output_record_start = target_gene['start'] - 1 - args.buffer_size
            output_record_end = target_gene['end'] - 1 + args.buffer_size
            output_record = record[output_record_start:output_record_end]
            Bio.SeqIO.write(output_record, sys.stdout, "fasta")
            
if __name__ == '__main__':
    main()
