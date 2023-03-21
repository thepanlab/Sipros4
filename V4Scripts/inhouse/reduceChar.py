#!/usr/bin/python


## Import Python package modules
import sys, getopt, warnings, os, re
from datetime import datetime, date, time
from collections import namedtuple
from collections import defaultdict
import csv
import math
import numpy as np


## Parse options
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
        output_filename = inputFileNameRoot + ".reducechar.txt"
    return (input_filename, output_filename)


def handleFileName(inputFileName, outputFileName) :
    inputFile = open(inputFileName)
    outputFile= open(outputFileName,"w")

    for eachLine in inputFile :
        eachLine = eachLine.strip()
        if (eachLine == ""):
            continue
        if (eachLine.startswith("#")):
            outputFile.write(eachLine+"\n")
            continue
        eachLine_List = eachLine.split("\t")
        outputLine = ""
        for eachField in eachLine_List :
            if (outputLine != "") :
                outputLine = outputLine + "\t"
            if len(eachField) > iFieldLimit :
                eachField = eachField[:iFieldLimit] + "..."
            outputLine = outputLine + eachField
        outputFile.write(outputLine+"\n")

    inputFile.close()
    outputFile.close()


## +------+
## | Main |
## +------+
def main(argv=None):

    # try to get arguments and error handling
    if argv is None:
        argv = sys.argv
        # parse options
        (inputFileName, outputFileName) = parse_options(argv)
        handleFileName(inputFileName, outputFileName)


global iFieldLimit 
iFieldLimit = 1000


## If this program runs as standalone, then go to main.
if __name__ == "__main__":
    main()
