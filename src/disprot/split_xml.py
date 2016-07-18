#!/usr/bin/env python

import sys
import xml
from xml.dom.minidom import parse, parseString
import xml.etree.ElementTree

__author__ = 'jurn'

def save_as_fasta(protein_id, sequence):
    filename = protein_id + ".fasta"
    fd = open(filename, "wb")
    fd.write(">" + protein_id + "\n")
    # array of lines
    # print sequence
    def split_lines(sequence_text):
        max_length = 60
        result = []
        remaining = sequence_text
        while len(remaining) > max_length:
            result.append(remaining[0:60])
            remaining = remaining[60:]
        result.append(remaining)
        return result
    sequence_lines = split_lines(sequence)
    for line in sequence_lines:
        fd.write(line + "\n")
    fd.close()
    print "wrote " + filename

def main():
    print "split xml"
    xml_file = sys.argv[1]

    tree = xml.etree.ElementTree.parse(xml_file)
    root = tree.getroot()
    for xml_protein in root.findall('{http://disprot.org/data/version_6.02/disprot_v6.02}protein'):
        protein_id = xml_protein.get('id')
        sequence = xml_protein.find('{http://disprot.org/data/version_6.02/disprot_v6.02}general/{http://disprot.org/data/version_6.02/disprot_v6.02}sequence').text
        save_as_fasta(protein_id, sequence)
        # write out a fasta file.

if __name__ == "__main__":
    main()

