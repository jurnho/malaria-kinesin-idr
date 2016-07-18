#!/usr/bin/env python

import sys
import xml
from xml.dom.minidom import parse, parseString
import xml.etree.ElementTree
from sets import Set

__author__ = 'jurn'



def main():
    print "disprot proteins with PDB entries"
    xml_file = sys.argv[1]

    tree = xml.etree.ElementTree.parse(xml_file)
    root = tree.getroot()
    total_count = 0
    for xml_protein in root.findall('{http://disprot.org/data/version_6.02/disprot_v6.02}protein'):
        protein_id = xml_protein.get('id')
        print protein_id
        xml_pdbs_ids = xml_protein.findall('d:regions/d:region/d:pdbs/d:pdb/d:id',
                                           namespaces={'d': 'http://disprot.org/data/version_6.02/disprot_v6.02'})
        if len(xml_pdbs_ids) > 0:
            total_count = total_count +1
        pdb_ids = Set()
        for x in xml_pdbs_ids:
            pdb_ids.add(x.text)
        print "  PDB: " + ",".join(pdb_ids)
    print "count of proteins with PDB entries: " + str(total_count)

if __name__ == "__main__":
    main()

