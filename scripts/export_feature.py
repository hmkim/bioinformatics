#!/usr/bin/env python

# Copyright 2017 by Junpeng Fan. All rights reserved.
#
# Last update 2017.4.25


import sys
import getopt

from os import path
from Bio import SeqIO
from Bio.SeqFeature import BeforePosition


Q_VALUE = {
    "gene": "gene",
    "coding sequences": "CDS",
    "ribosomal RNA": "rRNA",

    "5S ribosomal RNA": "5S",
    "5.8S ribosomal RNA": "5.8S",
    "16S ribosomal RNA": "16S",
    "18S ribosomal RNA": "18S",
    "23S ribosomal RNA": "23S",
    "28S ribosomal RNA": "28S",
    "internal transcribed spacer 1": "ITS1",
    "internal transcribed spacer 2": "ITS2"
}


def usage():
    q = ""
    for (key, value) in Q_VALUE.items():
        s = "    %-10s %s\n" % (value, key)
        q += s

    print("""usage: export_feature.py [options]

This script exports accurate seq from GenBank file according to information gaven.

Example: python export_feature.py -i Rhizopus_ITS.gb -t ITS1,5.8S,ITS2 -o ITS1-5.8S-ITS2.fa

Options:
-h --help    Print out this help message
-i --in      filename of GenBank file
-t --type    list of type included below. Multiple types should be separated by ","
    type      long name
%s-o --out     file name of output
    """ % q)


def read_seq(file):
    """
    read seq information from genbank file
    :param file: file name
    :return: seq object
    """

    # open file
    if path.exists(file):
        records = list(SeqIO.parse(file, 'gb'))
        return records
    else:
        raise ValueError("%s not exist" % file)


def extract(record, seq_type = "nucl", ftype = [], qual = [], qvalue = []):
    """
    extract seq or message from .gb
    :param record: Biopython Seq object
    :param seq_type: nucleotide or protein
    :param ftype:  feature type
    :param qual: qualifier name
    :param qvalue: qualifier value
    :return: fasta seq
    """

    print("Extract %s seq from record %s" % (", ".join(qvalue), record.id))
    r = ""

    for feature in record.features:
        # check featrue.type
        if ftype:
            if feature.type in ftype:
                pass
            else:
                continue
        else:
            pass

        start = feature.location.start
        end = feature.location.end
        strand = str(feature.strand)
        seq = record.seq[int(start):int(end)]

        # check strand
        if strand == "-1":
            seq = seq.reverse_complement()
        else:
            pass

        # check nucleotide or protein
        if seq_type == "prot":
            seq = seq.translate(feature.qualifiers.get("transl_table")[0])
        else:
            pass

        # check feature.qualifiers and qualifier value
        if qual:
            name = ""
            status = 0

            value_long = []  # long name of qvalue

            for n in qvalue:
                for key, value in Q_VALUE.items():
                    if n.lower() == value.lower():
                        value_long.append(key)
                    else:
                        pass
            #print(value_long)

            # get qual name and check status
            for q in qual:
                values = feature.qualifiers.get(q)

                # check qvalue
                if qvalue:
                    pass
                else:
                    value_long = values

                if values:
                    for value in values:
                        if value in value_long:
                            if value in Q_VALUE.keys():
                                name += "%s:%s " % (q, Q_VALUE[value])
                            else:
                                name += "%s:%s " % (q, value)
                            status = 1
                        else:
                            pass
                else:
                    pass

            if status == 1:
                # start position + 1
                if isinstance(start, BeforePosition):
                    start = "<" + str(int(start) + 1)
                else:
                    start += 1
                    pass

                print("%s %s-%s %s" % (name, start, end, strand))
                r += ">%s %s %s-%s %s\n%s\n" % (record.id, name, start, end, strand, seq)
            else:
                pass
        else:
            pass

    print("Record %s completed" % record.id)
    return r


def main():
    try:
        options, args = getopt.getopt(sys.argv[1:], "hi:t:o:", ["help", "in=", "type=", "out="])
    except getopt.error as x:
        print(x)
        sys.exit(0)

    #print(args)

    if len(args): # check extraneous arguments
        usage()
        sys.exit(0)
    else:
        pass

    in_name = ""
    type_list = ()
    out_name = ""

    for option, arg in options:
        if option in ("-h", "--help"):
            usage()
            sys.exit(0)
        elif option in ("-i", "--in"):
            in_name = arg
        elif option in ("-t", "--type"):
            type_list = arg.split(",")
        elif option in ("-o", "--out"):
            out_name = arg
        else:
            pass

    if in_name and type_list and out_name:
        output_file = open(out_name, 'w')
        records = read_seq(in_name)
        seq_type = "nucl"
        ftype = []
        qual = None
        qvalue = []

        if "CDS" in type_list:
            seq_type = "prot"
            ftype = ["CDS"]
            qual = "locus_tag"
        elif "gene" in type_list:
            ftype = ["gene"]
            qual = ["gene","locus_tag"]
        else:
            qual = ["product"]
            qvalue = type_list

        for record in records:
            seq = extract(record, seq_type, ftype, qual, qvalue)
            output_file.write(seq)
    else:
        raise ValueError("Options error. ")


if __name__ == '__main__':
    main()
