import getopt, sys
from urllib import urlencode
import cookielib, urllib2, os, re, copy, string, operator


def NewSearchName(currentLine, fileId) :
    allInfo = currentLine.split("=")
    sSearchName = allInfo[1]
    sSearchName = sSearchName.strip()
    if ((sSearchName == "Null") or (sSearchName == "null") or (sSearchName == "NULL")) :
        sSearchName = "subdb"+str(fileId)
    else :
        sSearchName = sSearchName + "_subdb"+str(fileId)
    return "Search_Name = "+sSearchName


def NewDB(currentLine, fileId, output_dir = "") :
    allInfo = currentLine.split("=")
    filePath = allInfo[1]
    filePath = filePath.strip()
    (pathRoot, pathExt) = os.path.splitext(filePath)
    if output_dir == "":
        return  "FASTA_Database = "+pathRoot+"_subfile"+str(fileId)+".fasta"
    else:
        # drive, path_and_file = os.path.splitdrive(pathRoot)
        # path, file = os.path.split(path_and_file)
        # return "FASTA_Database = "+os.path.join(output_dir, file+"_subfile"+str(fileId)+".fasta")
        return "FASTA_Database = "+output_dir+"_subfile"+str(fileId)+".fasta"


OriginalConfigureFileName = sys.argv[1]
OriginalConfigureFile = open(OriginalConfigureFileName)
outputpath            = sys.argv[2]
filenum               = int (sys.argv[3])
dboutputpath = ""
if len(sys.argv) == 5:
    dboutputpath = sys.argv[4]

ConfigureSubFiles = []

for i in range (filenum) :
        #os.popen("mkdir "+outputpath+str(i))
        (OriginalConfigureFileNameRoot, OriginalConfigureFileNameExt) = os.path.splitext(OriginalConfigureFileName)
        curruentConfigureSubFile = open(outputpath+os.sep+os.path.basename(OriginalConfigureFileNameRoot)+"_subdb"+str(i)+".cfg", "w")
        ConfigureSubFiles.append(curruentConfigureSubFile)

OriginalConfigureContent = OriginalConfigureFile.readlines()

for i in range (filenum) :
    for eachLine in OriginalConfigureContent :
        currentLine = eachLine.strip()
        if currentLine.startswith("Search_Name") :
            newLine = NewSearchName(currentLine, i)
            ConfigureSubFiles[i].write(newLine+"\n")
        elif currentLine.startswith("FASTA_Database") :
            newLine = NewDB(currentLine, i, dboutputpath)
            ConfigureSubFiles[i].write(newLine+"\n")
        else :
            ConfigureSubFiles[i].write(currentLine+"\n")

# close files
for i in range (filenum) :
        ConfigureSubFiles[i].close()

OriginalConfigureFile.close()