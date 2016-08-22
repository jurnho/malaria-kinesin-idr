#!/usr/bin/env python

# read the plasmo svm predictions, extract the results of a protein, and convert that into a Jalview annotation so
# we can visualise the result.

import sys
import csv

__author__ = 'jurn'


# convert a list of sorted positions into a list of ranges
def to_ranges(disordered_residue_positions):
    ranges = []
    current_range = None
    for p in disordered_residue_positions:
        if current_range is None:
            current_range = {
                'from': p,
                'to': p
            }
            continue
        elif current_range['to'] == p - 1:
            # extend the current range +1.
            current_range['to'] = p
        else:
            # start next
            ranges.append(current_range)
            current_range = {
                'from': p,
                'to': p
            }
    # add last
    if current_range is not None:
        ranges.append(current_range)
    return ranges


def main():
    if len(sys.argv) != 3:
        print "usage: " + sys.argv[0] + " [plasmodb_svm_results.csv] [motor_domains.csv]"
        return
    csv_results_file = sys.argv[1]
    motor_domains_file = sys.argv[2]

    # state held in here
    annotations = {}
    # add disorder predictions
    with open(csv_results_file, 'rb') as csv_file:
        csv_reader = csv.DictReader(csv_file)
        for row in csv_reader:
            # assume that position is already ordered.
            protein_id = row['protein_id']
            if protein_id not in annotations:
                annotations[protein_id] = {
                    "previous_position": 0,
                    "disordered_residue_positions": [],
                    "distances": []
                }
            protein_annotation = annotations[protein_id]

            distance = row['plasmo_svm_decision_values']
            disorder = row['plasmo_svm.disorder.smoothed']
            position = int(row['position'])
            if position != protein_annotation['previous_position'] + 1:
                raise Exception("missing residue at position: " + str(protein_annotation['previous_position']))
            protein_annotation['previous_position'] = int(position)
            protein_annotation['distances'].append(distance)
            if disorder == 'Y':
                protein_annotation["disordered_residue_positions"].append(position)
    # add motor domain info
    with open(motor_domains_file, 'rb') as csv_file:
        csv_reader = csv.DictReader(csv_file)
        for row in csv_reader:
            protein_id = row['protein id']
            motor_domain_start = row['motor domain start']
            motor_domain_end = row['motor domain end']
            protein_annotation = annotations[protein_id]
            protein_annotation['motor_domain_start'] = motor_domain_start
            protein_annotation['motor_domain_end'] = motor_domain_end
    # dump annotations
    for protein_id, prediction_info in annotations.iteritems():
        print protein_id
        output_file_name = protein_id + ".annotation"
        with open(output_file_name, 'w') as annotation_file:
            annotation_file.write("JALVIEW_ANNOTATION\n")
            annotation_file.write("SEQUENCE_REF\t" + protein_id + "\n")
            disordered_residue_ranges = to_ranges(prediction_info["disordered_residue_positions"])
            group_counter = 1
            for disordered_residue_range in disordered_residue_ranges:
                group_name = "Region_" + str(group_counter)
                group_counter += 1
                _from = str(disordered_residue_range['from'])
                to = str(disordered_residue_range['to'])
                annotation_file.write("SEQUENCE_GROUP\t" + group_name + "\t" + _from + "\t" +
                                      to + "\t*\n")
                annotation_file.write(
                        "PROPERTIES\t" + group_name + "\tdescription=disordered region\t"
                                                      "outlineColour=blue\tcolour=yellow\n")
            annotation_file.write("LINE_GRAPH\tSVM decision value\t" + "|".join(prediction_info["distances"]) + "\n")
            annotation_file.write("COLOUR\tSVM decision value\tblue\n")
            if 'motor_domain_start' in prediction_info:
                annotation_file.write('NO_GRAPH\tMotor Domain\t')
                start = prediction_info['motor_domain_start']
                # gaps, then H*
                annotation_file.write("|" * (int(start) -1))
                # move the label a bit further on.. jalview doesn't align it properlyy otherwise.
                annotation_file.write("H|H|H|H|H|H,motor domain|")
                end = prediction_info['motor_domain_end']
                motor_domain_length = int(end) - int(start)
                annotation_file.write("H|" * (motor_domain_length-6))

if __name__ == "__main__":
    main()
