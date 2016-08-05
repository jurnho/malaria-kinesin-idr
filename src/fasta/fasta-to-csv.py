#!/usr/bin/env python

import sys
import re
import csv

filename = sys.argv[1]
f = open (filename, 'r')
output_filename = re.sub('.fasta', '.csv', filename)

protein_id = ''
description = ''
sequence = ''
for line in f:
    line = line.rstrip()
    if line.startswith('>'):
        # print "header: " + line
        fields = line.split(' | ')
        protein_id = fields[0].replace('>','')
        # print fields
        description = fields[2]
    else:
        sequence = sequence + line

# print protein_id
# print sequence
sequence_length = len(sequence)
# print sequence_length
# print output_filename
with open(output_filename, 'wb') as csvfile:
    csvwriter = csv.writer(csvfile)
    csvwriter.writerow(["protein_id", "position", "residue"])
    for i in range(0, sequence_length):
        # print i
        residue = sequence[i]
        position = i+1
        csvwriter.writerow([protein_id, position, residue])