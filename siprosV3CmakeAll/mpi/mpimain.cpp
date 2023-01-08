#include <mpi.h>
#include <iostream>
#include <string>
#include <cstdlib>


#include "directoryStructure.h"
#include "proNovoConfig.h"
#include "ms2scanvector.h"

#define WORKTAG    1
#define DIETAG     2


using namespace std;

struct unit_of_workload_t
{
    string sFT2Filename;
    string sConfigureFilename;
    string sOutputDirectory;
};


void searchFT2Files (vector<string> & vsFT2Filenames, const string & sWorkingDirectory, bool bScreenOutput)
{
    int i, iFileNum;
    DirectoryStructure working_dir(sWorkingDirectory);
    working_dir.setPattern(".ft2");
    working_dir.getFiles(vsFT2Filenames);
    working_dir.setPattern(".FT2");
    working_dir.getFiles(vsFT2Filenames);
    working_dir.setPattern(".ms2");
    working_dir.getFiles(vsFT2Filenames);
    working_dir.setPattern(".MS2");
    working_dir.getFiles(vsFT2Filenames);
    iFileNum = (int) vsFT2Filenames.size();
    if (iFileNum == 0)
    {
	cerr<<"no scan file in the working directory"<<endl;
	exit(1);
    }
    for (i=0; i< iFileNum; i++)
	vsFT2Filenames.at(i) = sWorkingDirectory + ProNovoConfig::getSeparator() + vsFT2Filenames.at(i);
}

void searchConfigureFiles (vector<string> & vsConfigureFilenames, const string & sConfigFileDirectory, bool bScreenOutput)
{
    int i, iFileNum;
    DirectoryStructure working_dir(sConfigFileDirectory);
    working_dir.setPattern(".cfg");
    working_dir.getFiles(vsConfigureFilenames);
    working_dir.setPattern(".CFG");
    working_dir.getFiles(vsConfigureFilenames);

    iFileNum = (int) vsConfigureFilenames.size();
    if (iFileNum == 0)
    {
	cerr<<"no configure file in the directory"<<endl;
	exit(1);
    }
    for (i=0; i< iFileNum; i++)
	vsConfigureFilenames.at(i) = sConfigFileDirectory + ProNovoConfig::getSeparator() + vsConfigureFilenames.at(i);
}

/* 
 * Parse command line arguments
 * Populate vsFT2Filenames
 * Set up SiprosConfig
 */

void initializeArguments(int argc, char **argv, vector<string> & vsFT2Filenames, 
			 string & sWorkingDirectory, vector<string> & vsConfigureFilenames, 
			 string & sSingleWorkingFile, string & sOutputDirectory,
			 bool & bScreenOutput)
// Under MPI mode, a user can specify configure file directory by specifying -g
// If no configre file or directory is specified, configure file directory is working directory 
{
    int i;

    string sConfigFileDirectory, sConfigFilename;
    // Grab command line arguments
    vector<string> vsArguments;
    
    sWorkingDirectory  = "";
    sConfigFilename    = "";
    sSingleWorkingFile = "";
    sOutputDirectory   = "";
    sConfigFileDirectory = "";
    bScreenOutput      = true;    
    
    while(argc--) 
	vsArguments.push_back(*argv++);
    for(i = 1; i <= (int)vsArguments.size()-1; i++)
	if(vsArguments[i] == "-w") 
		sWorkingDirectory = vsArguments[++i]; 
	else if (vsArguments[i] == "-c")
		sConfigFilename = vsArguments[++i];
	else if (vsArguments[i] == "-f")
		sSingleWorkingFile = vsArguments[++i];
	else if (vsArguments[i] == "-o")
		sOutputDirectory = vsArguments[++i];
	else if (vsArguments[i] == "-g")
		sConfigFileDirectory = vsArguments[++i];
	else if (vsArguments[i] == "-s")
		bScreenOutput = false;
	else if ((vsArguments[i] == "-h") || (vsArguments[i] == "--help"))
	{
		cout << "Usage: -w WorkingDirectory -c ConfigurationFile, -f: A single MS2 or FT2 file to be processed" << endl;
		cout << "If configuration file is not specified, Sipros will look for SiprosConfig.cfg in the directory of FT2 files" << endl;
		cout << "-o output directory. If not specified, it is the same as that of the input scan file,"<<endl;
		cout << "-g configure file directory. -s silence all standard output."<<endl;
		exit(0);
	}else 
	{
		cerr << "Unknown option " << vsArguments[i] << endl << endl; 
		exit(1);
	}
    if ((sWorkingDirectory == "") && (sSingleWorkingFile == "") )
	sWorkingDirectory = ".";
	
    if ((sWorkingDirectory != "") && (sSingleWorkingFile != ""))
    {
	cerr << "Either a input scan file or the directory of input scan files needs to be specified" <<endl;
	exit(1);
    }
    if ((sConfigFilename == "") && (sConfigFileDirectory == ""))
	//sConfigFilename = sWorkingDirectory + ProNovoConfig::getSeparator() + "SiprosConfig.cfg";
	// Without specifying configure file and configure file directory, the default configure file directory
	// is working directory. In the openmp only version, the default configure file is SiprosConfig.cfg
	sConfigFileDirectory = sWorkingDirectory;

    if (sConfigFileDirectory != "")
        searchConfigureFiles (vsConfigureFilenames, sConfigFileDirectory, bScreenOutput);
    else
        vsConfigureFilenames.push_back(sConfigFilename);

    if (sSingleWorkingFile != "")
	vsFT2Filenames.push_back(sSingleWorkingFile);
    else
	searchFT2Files (vsFT2Filenames, sWorkingDirectory, bScreenOutput);
    if ((sOutputDirectory == "") && (sWorkingDirectory != ""))
	sOutputDirectory = sWorkingDirectory;
    

}


void handleScan(const string & sFT2filename, const string & sOutputDirectory, const string & sConfigFilename, bool bScreenOutput)
{
    MS2ScanVector * pMainMS2ScanVector = new MS2ScanVector(sFT2filename, sOutputDirectory, sConfigFilename, bScreenOutput);
    
    if (bScreenOutput)
	cout<<"Reading MS2 scan file "<<sFT2filename<<endl;
    
    if(!pMainMS2ScanVector->loadFT2file())
	cerr << "Error: Failed to load file: " << sFT2filename << endl;
    else
    {
    // search all MS2 scans and write output to a file
	pMainMS2ScanVector->startProcessing();
    }
    delete pMainMS2ScanVector;//free memory of vpAllMS2Scans   
}


void MasterProcess(const vector <unit_of_workload_t> & vWorkload, bool bScreenOutput)
{
    size_t i, workloadSize, iBounderOfProcess;
    int  currentWorkId; //unit id of vWorkLoad
    int  iNumberOfProcessors, iNumberOfSlaves;
    int  result;
    MPI_Status status;
    MPI_Comm_size(MPI_COMM_WORLD, &iNumberOfProcessors);     /* get number of processes */
    workloadSize = vWorkload.size();

    iNumberOfSlaves = iNumberOfProcessors -1;
    iBounderOfProcess = ((workloadSize <= (size_t)iNumberOfSlaves) ? workloadSize : (size_t)iNumberOfSlaves) ;
    for (i=1; i<=iBounderOfProcess; i++)
    {
        currentWorkId = i-1;
        MPI_Send(&currentWorkId, 1, MPI_INT, i, WORKTAG, MPI_COMM_WORLD);
    }
    if ((int)workloadSize >  iNumberOfSlaves)
    {
        currentWorkId = iNumberOfSlaves;
        while (currentWorkId < (int) workloadSize)
        {
            MPI_Recv(&result, 1, MPI_INT, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
            MPI_Send(&currentWorkId, 1, MPI_INT, status.MPI_SOURCE, WORKTAG, MPI_COMM_WORLD);
            currentWorkId ++;
        }
    }
    /* Tell all the slaves to exit by sending an empty message with the DIETAG. */
    for (i=1; i<= (size_t)iNumberOfSlaves; i++)
        MPI_Send(0, 0, MPI_INT, i, DIETAG, MPI_COMM_WORLD);
    if (bScreenOutput)
	cout<<"Master process is done."<<endl;
}




void SlaveProcess(const vector <unit_of_workload_t> & vWorkload, bool bScreenOutput)
{
    MPI_Status status;
    int currentWorkId, myid;
    unit_of_workload_t currentWork;
    MPI_Comm_rank(MPI_COMM_WORLD, &myid);
    while(true)
    {
        MPI_Recv(&currentWorkId, 1, MPI_INT, 0, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
        if (status.MPI_TAG == DIETAG)
            break;
        currentWork = vWorkload.at(currentWorkId);
        //cout<<"slave id:"<<myid<<" sFT2name: "<<currentWork.sFT2Filename<<endl;
        // Load config file
        if(!ProNovoConfig::setFilename(currentWork.sConfigureFilename))
        {
	    cerr << "Could not load config file " << currentWork.sConfigureFilename << endl;
            exit(1);
        }
        handleScan(currentWork.sFT2Filename, currentWork.sOutputDirectory, currentWork.sConfigureFilename, bScreenOutput);
	if (bScreenOutput)
	    cout<<currentWork.sFT2Filename<<" and "<<currentWork.sConfigureFilename<<" is done by Slave process "<<myid<<endl;
        MPI_Send(0, 0, MPI_INT, 0, 0, MPI_COMM_WORLD);
    }
    if (bScreenOutput)
	cout<<"Slave process "<<myid<<" is done."<<endl;
}




int main(int argc, char **argv) 
{
    // A list of FT2/MS2 files to be searched
    vector<string> vsFT2Filenames;
    // A list of configure files  
    vector<string> vsConfigureFilenames;
    int i, j, myid;  
    string sWorkingDirectory, sConfigFilename, sSingleWorkingFile, sOutputDirectory, sConfigFileDirectory;
    bool bScreenOutput;
    unit_of_workload_t current_work;
    vector <unit_of_workload_t> vWorkload;
    MPI_Init(&argc,&argv);              /* starts MPI */
    initializeArguments(argc, argv, vsFT2Filenames, sWorkingDirectory, vsConfigureFilenames, sSingleWorkingFile,
			 sOutputDirectory, bScreenOutput);

    for (j = 0; j < (int) vsConfigureFilenames.size(); j++)
    {
        sConfigFilename = vsConfigureFilenames[j] ;
        
        // Process one FT2 file and one configure file  at a time  
        for(i = 0; i < (int) vsFT2Filenames.size(); i++)
        {
            current_work.sFT2Filename = vsFT2Filenames[i];
            current_work.sConfigureFilename = sConfigFilename;
            current_work.sOutputDirectory = sOutputDirectory;
            vWorkload.push_back(current_work);
	    //handleScan(vsFT2Filenames[i], sOutputDirectory, sConfigFilename);
        }
    }

    MPI_Comm_rank(MPI_COMM_WORLD, &myid);  /* get current process id */
    if (myid == 0 )
        MasterProcess(vWorkload, bScreenOutput);
    else
        SlaveProcess(vWorkload, bScreenOutput);
    MPI_Finalize();          /* let MPI finish up ... */
    //std::cout << "Hello, world!" << std::endl;
    return 0;
}
