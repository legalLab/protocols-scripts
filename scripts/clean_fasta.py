#! /usr/bin/env python3

#######################################
# A simple script to remove short sequences from a FASTA file
#
# FASTA cleanup
#
# Author:  Tomas Hrbek
# Email: hrbek@evoamazon.net
# Date: 26.06.2015
# Version: 1.0
#######################################

import argparse
from Bio import SeqIO

__author__ = 'legal'

 
parser = argparse.ArgumentParser(description='Script to remove short sequences from a fast file.')
parser.add_argument('input', help='input file name')
parser.add_argument('output', help='output file name')
parser.add_argument('cut_off', help='minimum fragment length, default 250 bp', type=int, nargs='?', default=250)
args = parser.parse_args()
    
n_reads = 0 #Initialize rad counter
sequences = [] #Setup an empty list of sequence lengths

input_handle = open(args.input, "rU")
output_handle = open(args.output, "w")

reads_iterator = SeqIO.parse(input_handle, "fasta")
seqs_iterator = (reads for reads in reads_iterator if len(reads) >= args.cut_off)

SeqIO.write(seqs_iterator, output_handle, "fasta")
output_handle.close()


## show values ##
print("Input file: {}".format(args.input))
print("Output file: {}".format(args.output))
print("Minimum read length: {}".format(args.cut_off))
