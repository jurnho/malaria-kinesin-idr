#!/usr/bin/env python

import os
import re
import sys
import csv

__author__ = 'jurn'


# converts stdout of a DisEMBL run into a CSV.
def save_as_csv(filename):
    if not filename.endswith(".txt"):
        raise Exception("expected file ending in .txt")
    output_filename = re.sub('.txt$', '_disembl.csv', filename)
    basename = os.path.basename(filename)
    protein_id = re.sub('.txt$', '', basename)
    print protein_id
    position = 0
    with open(output_filename, 'wb') as csvfile:

        csvwriter = csv.writer(csvfile)
        csvwriter.writerow(["protein_id", "predictor", "position", "score"])
        for line in open(filename):
            line = line.rstrip()
            if re.match("^[A-Z]", line):
                position = position + 1
                # print line
                fields = line.split("\t")
                # print fields
                csvwriter.writerow([protein_id, "DISEMBL_COILS", position, fields[1]])
                csvwriter.writerow([protein_id, "DISEMBL_REM465", position, fields[2]])
                csvwriter.writerow([protein_id, "DISEMBL_HOTLOOPS", position, fields[3]])


def main():
    # print "to csv"
    filename = sys.argv[1]

    save_as_csv(filename)


if __name__ == "__main__":
    main()
