#!/usr/bin/env python

import sys
import xml
from xml.dom.minidom import parse, parseString
import xml.etree.ElementTree
import json
import io
import re
import csv
import os
__author__ = 'jurn'

def save_as_csv(filename):
    basename = os.path.basename(filename)
    protein_id = re.sub(".diso", "", basename)
    output_filename = re.sub(".diso", "_" + "disopred" + ".csv", filename)
    with open(output_filename, 'wb') as csvfile:
        csvwriter = csv.writer(csvfile)
        csvwriter.writerow(["protein_id", "predictor", "position", "score"])
        position = 0
        for line in open(filename):
            line = line.rstrip().lstrip()
            if not line.startswith("#"):
                # print line
                position = position + 1
                fields = re.split("\s+", line)
                # print fields
                csvwriter.writerow([protein_id, "DISOPRED", position, fields[3]])


def main():
    # print "to json"
    filename = sys.argv[1]

    save_as_csv(filename)

if __name__ == "__main__":
    main()

