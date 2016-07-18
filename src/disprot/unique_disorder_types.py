#!/usr/bin/env python

import sys
import xml
from xml.dom.minidom import parse, parseString
import xml.etree.ElementTree

__author__ = 'jurn'

def main():
    print "unique disorder types"
    xml_file = sys.argv[1]

    tree = xml.etree.ElementTree.parse(xml_file)
    root = tree.getroot()
    for xml_protein in root.findall('{http://disprot.org/data/version_6.02/disprot_v6.02}protein'):
        protein_id = xml_protein.get('id')
        # sequence = xml_protein.find('{http://disprot.org/data/version_6.02/disprot_v6.02}general/{http://disprot.org/data/version_6.02/disprot_v6.02}sequence').text
        # save_as_fasta(protein_id, sequence)
        # write out a fasta file.
        print protein_id
        xml_regions = xml_protein.findall('{http://disprot.org/data/version_6.02/disprot_v6.02}regions/{http://disprot.org/data/version_6.02/disprot_v6.02}region')

        for xml_region in xml_regions:
            region_id = xml_region.get('id')
            type = xml_region.find('{http://disprot.org/data/version_6.02/disprot_v6.02}type').text
            print "  " + region_id + ":" + type

if __name__ == "__main__":
    main()

