# -*- coding: utf-8 -*-


import sys, getopt, warnings, os, re
from datetime import datetime, date, time


def parse_options(argv):
    
    opts, args = getopt.getopt(argv[1:], "hw:n:o:")

    working_dir = ""
    filenumber  = ""
    outputpath  = ""

    # Basic options
    for option, value in opts:
        if option in ("-h"):
            print("-w workingdirectory -n num_split -o out_directory")
            sys.exit(1)
        if option in ("-w"):
            working_dir = value
        if option in ("-n",):   
            filenumber = value
        if option in ("-o"):    
            outputpath = value

    if filenumber == "" or outputpath == "" or working_dir == "" :
        print("please specify -w workingdirectory -n num_split -o out_directory")
        sys.exit(1)
    
    FT2_filename_list = get_file_list_with_ext(working_dir, ".ft2")
    MS2_filename_list = get_file_list_with_ext(working_dir, ".ms2")
    
    Scans_filename_list = FT2_filename_list + MS2_filename_list 

    return [Scans_filename_list, filenumber, outputpath]

## Get file(s) list in working dir with specific file extension
def get_file_list_with_ext(working_dir, file_ext):

    # define sipros file extension 
    file_list = []

    # working directory
    if os.path.exists(working_dir):
        for file_name in os.listdir(working_dir):

            # check the file extension
            if file_name.lower().endswith(file_ext):
                file_path_name = os.path.join(working_dir, file_name)
                file_list.append(file_path_name)

        file_list = sorted(file_list)

    else:
        psys.stderr.write("\nCannot open working directory {}".format(working_dir))
        sys.exit(1)

    return file_list

## +------+
## | Main |
## +------+
def main(argv=None):

    # try to get arguments and error handling
    if argv is None:
		argv = sys.argv
    
    # parse options
    [Scans_filename_list, filenumber, outputpath] = parse_options(argv)
    
    for eachFileName in Scans_filename_list :
        os.system("python separate_scans.py -f " + eachFileName + " -o " + outputpath + " -n " +filenumber)

## If this program runs as standalone, then exit.
if __name__ == "__main__":
    sys.exit(main())




