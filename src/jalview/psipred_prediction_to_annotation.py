#!/usr/bin/env python

import sys
__author__ = 'jurn'


# convert a psipred 4.0 ss2 secondary structure prediction into a jalview annotation
def main():
    if len(sys.argv) != 3:
        print "usage: " + sys.argv[0] + " [psipred ss2 file] [jalview annotation file]"
        return
    ss2_prediction_file = sys.argv[1]
    output_file_name = sys.argv[2]
    with open(output_file_name, 'w') as annotation_file:
        annotation_file.write("JALVIEW_ANNOTATION\n")
        annotation_file.write("NO_GRAPH\tPSIPRED\t")

        f = open(ss2_prediction_file, 'r')
        for line in f:
            line = line.strip()
            if line.startswith("#"):
                continue
            fields = line.split(" ")
            if len(fields) < 3:
                continue
            sse = fields[2]
            if sse == 'C':
                annotation_file.write("|")
            elif sse == 'H':
                annotation_file.write("H|")
            elif sse == 'E':
                annotation_file.write("E|")
            else:
                raise Exception("unhandled: " + sse)

        annotation_file.write("\n")

if __name__ == "__main__":
    main()
