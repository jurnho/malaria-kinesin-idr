#!/usr/bin/env python

import sys
import xml
from xml.dom.minidom import parse, parseString
import xml.etree.ElementTree
import json
import io

__author__ = 'jurn'

# http://www.disprot.org/disorderChar.php
def save_as_json(protein_id, sequence, xml_protein):
    filename = protein_id + ".json"
    data = {
        "id": protein_id,
        "sequence": sequence,
        "regions" : [],
        "predictor": "disprot"
    }

    for region in xml_protein.findall('{http://disprot.org/data/version_6.02/disprot_v6.02}regions/{http://disprot.org/data/version_6.02/disprot_v6.02}region'):
        region_type = region.find('{http://disprot.org/data/version_6.02/disprot_v6.02}type').text
        region_start = region.find('{http://disprot.org/data/version_6.02/disprot_v6.02}start').text
        region_end = region.find('{http://disprot.org/data/version_6.02/disprot_v6.02}end').text
        data["regions"].append({
            "type": region_type,
            "start": region_start,
            "end": region_end
        })
        print region_type
    fd = io.open(filename, "wb")
    json.dump(data, fd)
    fd.close()


def main():
    print "split xml"
    xml_file = sys.argv[1]

    tree = xml.etree.ElementTree.parse(xml_file)
    root = tree.getroot()
    for xml_protein in root.findall('{http://disprot.org/data/version_6.02/disprot_v6.02}protein'):
        protein_id = xml_protein.get('id')
        sequence = xml_protein.find('{http://disprot.org/data/version_6.02/disprot_v6.02}general/{http://disprot.org/data/version_6.02/disprot_v6.02}sequence').text
        save_as_json(protein_id, sequence, xml_protein)
        # write out a json file.

if __name__ == "__main__":
    main()

