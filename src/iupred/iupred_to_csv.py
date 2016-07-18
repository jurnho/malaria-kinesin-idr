#!/usr/bin/env python

import sys
import csv
import xml
from xml.dom.minidom import parse, parseString
import xml.etree.ElementTree
import json
import io
import re
__author__ = 'jurn'

def save_as_csv(prediction_type, filename):
    sequence = ''
    regions = []

    # contains protein id
    last_comment = ''
    if not filename.endswith(".txt"):
        raise Exception("expected file ending in .txt")
    predictor = "iupred_" + prediction_type
    output_filename = re.sub('.txt$', '_' + predictor + '.csv', filename)
    with open(output_filename, 'wb') as csvfile:
        csvwriter = csv.writer(csvfile)
        csvwriter.writerow(["protein_id", "predictor", "position", "score"])

        for line in open(filename):
            # print "line: " + line
            line = line.rstrip()
            line = line.lstrip()
            if line.startswith("#"):
                last_comment = re.sub("# ", '', line)
                continue
            protein_id = last_comment
            fields = re.split('\\s+', line)
            # print fields
            sequence += fields[1]
            score = float(fields[2])
            position = fields[0]
            csvwriter.writerow([protein_id, predictor, position, score])


def main():
    if len(sys.argv) != 3:
        print "usage: " + sys.argv[0] + " [long|short] iupred_output_file"
    prediction_type = sys.argv[1]
    filename = sys.argv[2]

    save_as_csv(prediction_type, filename)

if __name__ == "__main__":
    main()

