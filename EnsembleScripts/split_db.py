import getopt, sys
from urllib import urlencode
import cookielib, urllib2, os, re, copy, string, operator

OriginalFastaFileName = sys.argv[1]
outputpath            = sys.argv[2]
filenum               = int (sys.argv[3])

FastaSubFiles = []

for i in range (filenum) :
        #os.popen("mkdir "+outputpath+str(i))
        (OriginalFastaFileNameRoot, OriginalFastaFileNameExt) = os.path.splitext(OriginalFastaFileName)
        curruentFastaSubFile = open(outputpath+os.sep+os.path.basename(OriginalFastaFileNameRoot)+"_subfile"+str(i)+".fasta", "w")
        FastaSubFiles.append(curruentFastaSubFile)

seqCount = -1

OriginalFastaFile = open(OriginalFastaFileName, "r")
while 1:
	sLine = OriginalFastaFile.readline()
	if not sLine:
		break
	if sLine[0] == ">":
		seqCount = seqCount + 1
		fileId = seqCount % filenum
		FastaSubFiles[fileId].write(sLine)
	else:
		fileId = seqCount % filenum
		FastaSubFiles[fileId].write(sLine)

# close files
OriginalFastaFile.close()
for i in range (filenum) :
        FastaSubFiles[i].close()
        
drive, path_and_file = os.path.splitdrive(OriginalFastaFileName)
path, file = os.path.split(path_and_file)
(pathRoot, pathExt) = os.path.splitext(file)
print os.path.join(outputpath, pathRoot)