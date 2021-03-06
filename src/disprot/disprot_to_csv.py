#!/usr/bin/env python

import csv
import os
import sys
import xml
import xml.etree.ElementTree
from collections import namedtuple
from ftplib import FTP

__author__ = 'jurn'

namespaces = {'d': 'http://disprot.org/data/version_6.02/disprot_v6.02'}

def is_disordered(position, regions):
    for r in regions:
        # print position
        # print r
        if r['start'] <= position <= r['end']:
            return {
                "disordered": 'Y',
                "type": r['type']
            }
    return {
        "disordered": 'N',
        "type": None
    }


PdbChain = namedtuple('PdbChain', ['id', 'chain'])


def get_pdb(xml_region):
    if (xml_region is None):
        return PdbChain("", "")
    pdb_id = xml_region.find('d:id', namespaces).text
    pdb_chain = xml_region.find('d:chain', namespaces).text
    return PdbChain(pdb_id, pdb_chain)


# secondary structure element, amino acid (residue).
DsspInfo = namedtuple('DsspInfo', ['sse', 'aa'])
get_dssp_see_cache = {}

def get_dssp_sse(pdb_id, chain, position, aa):
    # mem cache
    if pdb_id == "":
        return DsspInfo("", "")
    if pdb_id in get_dssp_see_cache:
        c = get_dssp_see_cache[pdb_id]
        # file not found
        if c is None:
            return DsspInfo("", "")
        else:
            pos_key = str(position)+"." + chain
            if pos_key in get_dssp_see_cache[pdb_id]:
                return get_dssp_see_cache[pdb_id][pos_key]
            else:
                return  DsspInfo("", "")

    # cachedir
    if not os.path.isdir('dssp'):
        os.mkdir('dssp')
    # check if already downloaded
    file_name = pdb_id.lower() + ".dssp"
    local_file_name = "dssp/" + file_name
    if not os.path.isfile(local_file_name):
        # download dssp via ftp
        if pdb_id == "":
            return DsspInfo("", "")
        url = "ftp://ftp.cmbi.ru.nl/pub/molbio/data/dssp/" + pdb_id.lower() + ".dssp"
        print url
        ftp = FTP('ftp.cmbi.ru.nl')
        ftp.set_debuglevel(2)
        ftp.login()
        ftp.cwd('/pub/molbio/data/dssp/')
        retr_command = 'RETR ' + file_name
        try:
            ftp.retrbinary(retr_command, open(local_file_name, 'wb').write)
        except:
            # may be superseded, skip
            print file_name + " not found"
            open(local_file_name, 'wb').write("")
            ftp.quit()
            return DsspInfo("", "")
        ftp.quit()
    # parse it
    # open local file name
    data = open(local_file_name, "rb").read()
    if len(data) == 0:
        # set as empty in cache
        get_dssp_see_cache[pdb_id] = None
        return DsspInfo("", "")
    print "pdb_id: " + pdb_id
    print "chain: " + chain
    print "position: " + str(position)
    print "aa: " + aa
    # raise Exception("testing")
    dssp_dict = parse_dssp(data)
    get_dssp_see_cache[pdb_id] = dssp_dict

    pos_key = str(position) + "." + chain

    # print dssp_dict
    # print pdb_id
    if pos_key in get_dssp_see_cache[pdb_id]:
        return dssp_dict[pos_key]
    else:
        return DsspInfo("", "")

def parse_dssp(dssp_data):
    lines = dssp_data.split("\n")
    in_data_records = False
    result = {}
    for line in lines:
        # print line
        if "RESIDUE AA STRUCTURE" in line:
            in_data_records = True
            # print "jurN"
            # print line
            continue
        if in_data_records:
            # print line[:17]
            residue_number = line[6:10].lstrip()
            chain = line[11:12]
            aa = line[13:14]
            sse = line[16:17].lstrip()
            # print "x" + residue_number + "x"
            # print "x" + chain + "x"
            # print "x" + sse + "x"
            result[residue_number + "." + chain] = DsspInfo(sse, aa)
    return result


def main(options):
    xml_file = sys.argv[1]

    tree = xml.etree.ElementTree.parse(xml_file)
    root = tree.getroot()

    if not options["split"]:
        csvfile = open("disprot.csv", 'wb')
        csvwriter = csv.writer(csvfile)
        csvwriter.writerow(
                ["protein_id", "uniprot_id", "name", "position", "residue", "disordered", "type", "pdb id",
                 "pdb chain", "dssp ss", "dssp aa"])

    for xml_protein in root.findall('d:protein', namespaces):
        protein_id = xml_protein.get('id')
        if options["split"]:
            # start a new file
            csvfile = open(protein_id + "_disprot.csv", 'wb')
            csvwriter = csv.writer(csvfile)
            csvwriter.writerow(
                    ["protein_id", "uniprot_id", "name", "position", "residue", "disordered", "type", "pdb id",
                     "pdb chain", "dssp ss", "dssp aa"])
        sequence = xml_protein.find(
                'd:general/d:sequence', namespaces).text
        uniprot_id = xml_protein.find(
                'd:general/d:uniprot', namespaces).text
        name = xml_protein.find(
                'd:general/d:name', namespaces).text
        xml_pdbs = xml_protein.findall('d:regions/d:region/d:pdbs/d:pdb',namespaces)
        if len(xml_pdbs) > 0:
            first_pdb_xml = xml_pdbs[0]
        else:
            first_pdb_xml = None
        pdb_entry = get_pdb(first_pdb_xml)

        print protein_id
        xml_regions = xml_protein.findall('d:regions/d:region', namespaces)

        regions = []
        for xml_region in xml_regions:
            # region_id = xml_region.get('id')
            type = xml_region.find('d:type', namespaces).text
            if type == "Ordered":
                continue
            # print "  " + region_id + ":" + type
            start = xml_region.find('d:start', namespaces).text
            start = int(start)
            end = xml_region.find('d:end', namespaces).text
            end = int(end)
            regions.append({
                "start": start,
                "end": end,
                "type": type
            })
        position = 0
        for aa in sequence:
            position = position + 1
            # is ordered
            disordered_result = is_disordered(position, regions)
            pdb_id = pdb_entry.id
            pdb_chain = pdb_entry.chain
            dssp_result = get_dssp_sse(pdb_id, pdb_chain, position, aa)
            sse = dssp_result.sse
            dssp_aa = dssp_result.aa
            csvwriter.writerow([protein_id, uniprot_id, name, position, aa, disordered_result['disordered'],
                                disordered_result['type'],
                                pdb_id, pdb_chain, sse, dssp_aa])
            # print aa
        if options["split"]:
            csvfile.close()

if __name__ == "__main__":
    if len(sys.argv) == 1:
        print "Usage: " + sys.argv[0] + " [--split] DISPROT_XML_FILE"
        print
        print "convert the DisProt XML database into a CSV file named disprot.csv"
        print
        print "options:"
        print "  split    instead of writing a single disprot.csv as output, create one csv per protein."
        print "           Files will be created in the format <protein_id>_disprot.csv"
        exit(0)
    else:
        options = {
            "split": False
        }
        for arg in sys.argv[1:]:
            if arg == '--split':
                options['split'] = True
            else:
                options['disprot_xml_file'] = arg
        main(options)

# protein_id, position, residue, disordered Y/N, type
