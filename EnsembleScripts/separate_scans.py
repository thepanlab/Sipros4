 # -*- coding: utf-8 -*-

 
import getopt, sys 
from urllib import urlencode
import cookielib, urllib2, os, re, copy, string, operator

def divide_ms2_file(FT2FileName, output_dir_str, filenum):
		
	FT2file = open(FT2FileName)
	filenum = int(filenum) 
	
	(pathRoot, pathExt) = os.path.splitext(FT2FileName)
	FT2subfiles = []
	
	for i in range (filenum) :
		#os.popen("mkdir "+output_dir_str+str(i))
		(FT2FileNameRoot, FT2FileNameExt) = os.path.splitext(FT2FileName)
		curFT2subfile = open(os.path.join(output_dir_str, os.path.basename(FT2FileNameRoot) + "_subfile" + str(i) + pathExt), "w")
		FT2subfiles.append(curFT2subfile)
	
	scannum = 0
	fileid  = -1
	
	for eachline in FT2file :
		if (eachline.startswith("H")):
			continue
		if (eachline.startswith("S")):
			scannum = scannum + 1
			fileid  = scannum % filenum
		FT2subfiles[fileid].write(eachline)
	
	# close files	
	for i in range (filenum) :
		FT2subfiles[i].close()
		
	FT2file.close()


def parse_options(argv):
    
    opts, args = getopt.getopt(argv[1:], "hf:n:o:")

    ms2_file_str = ""
    num_split  = ""
    output_dir_str  = ""

    # Basic options
    for option, value in opts:
        if option in ("-h"):
            print("-f ms2_file -n num_split -o out_directory")
            sys.exit(1)
        if option in ("-f"):
            ms2_file_str = value
        if option in ("-n",):   
            num_split = value
        if option in ("-o"):    
            output_dir_str = value

    if num_split == "" or output_dir_str == "" or ms2_file_str == "" :
        print("please specify -f ms2_file -n num_split -o out_directory")
        sys.exit(1)

    return [ms2_file_str, num_split, output_dir_str]

## +------+
## | Main |
## +------+
def main(argv=None):
	
	if argv is None:
		argv = sys.argv
		
	print(' '.join(sys.argv))
	
	# parse options
	[ms2_file_str, num_split, output_dir_str] = parse_options(argv)
	
	# split file
	divide_ms2_file(ms2_file_str, output_dir_str, num_split)
	
	print('Done.')

## If this program runs as standalone, then exit.
if __name__ == "__main__":
    sys.exit(main())
