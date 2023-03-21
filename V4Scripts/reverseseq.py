#!/usr/bin/python


## Import Python package modules
import sys, getopt, warnings, os, re
from datetime import datetime, date, time
from collections import namedtuple
from collections import defaultdict
import csv
import math
import numpy as np
from Bio import SeqIO
from Bio.Seq import Seq

import parseconfig
import sipros_post_module


def parse_options(argv):

    opts, args = getopt.getopt(argv[1:], "hi:o:",
                                    ["help",
                             	     "input-file",
	                			     "output-file"])

    output_filename = ""
    input_filename  = ""

    # Basic options
    for option, value in opts:
        if option in ("-h", "--help"):
            print "-i input-file, -o output-file"
            sys.exit(1)
        if option in ("-i", "--input-file"):
            input_filename = value
        if option in ("-o", "--output-file"):
            output_filename = value

    if (input_filename == "") :
        print "Please specify -i"
        sys.exit(1)
    if (output_filename == "") :
        (inputFileNameRoot, inputFileNameExt) = os.path.splitext(input_filename)
        output_filename = inputFileNameRoot + "_CFR" + inputFileNameExt
    return (input_filename, output_filename)


def ReverseSeq(inputFileName, outputFileName) :
    outputFile = open(outputFileName, "w")
    for record in SeqIO.parse(inputFileName, "fasta" ) :
        currentSeq = str(record.seq)
        if not(currentSeq.isalpha()) :
            print "Remove non-alphabetic characters for sequence "+record.id
            currentSeq = re.sub(r'\W+', '', currentSeq)
        outputFile.write(">"+record.description+"\n")
        outputFile.write(currentSeq+"\n")
        outputFile.write(">Rev_"+record.id+" Decoy\n")
        outputFile.write(currentSeq[::-1]+"\n")


    outputFile.close()


def main(argv=None):

    # try to get arguments and error handling
    if argv is None:
        argv = sys.argv
        # parse options
        (inputFileName, outputFileName) = parse_options(argv)
        ReverseSeq(inputFileName, outputFileName)


## If this program runs as standalone, then go to main.
if __name__ == "__main__":
    main()




