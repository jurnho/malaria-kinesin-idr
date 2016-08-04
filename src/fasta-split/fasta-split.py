#!/usr/bin/env python

# split a single fasta file with multiple sequences into individual fasta files, using their id as part of the filename
from Bio import SeqIO
import sys

filename = sys.argv[1]
handle = open(filename, "rU")
for record in SeqIO.parse(handle, "fasta") :
    print record.id
    # write out per file
    out = open(record.id + ".fasta" , "w")
    SeqIO.write(record, out, "fasta")
    out.close()
handle.close()
