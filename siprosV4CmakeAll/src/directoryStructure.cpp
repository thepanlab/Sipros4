
#include "directoryStructure.h"
#include <sys/types.h>
#include <dirent.h>

DirectoryStructure::DirectoryStructure( const string & sDir )
{
	iFileCount = 0;
	sDirectory = sDir;
}

DirectoryStructure::DirectoryStructure()
{
	iFileCount = 0;
}

DirectoryStructure::~DirectoryStructure()
{

}

void DirectoryStructure::setDirectory( const string & sDir )
{
	sDirectory = sDir;
}

void DirectoryStructure::setPattern( string sPattern )
{

	vsFileList.clear();
	
	DIR *d;
	struct dirent *dir;
	d = opendir(sDirectory.c_str());
	if (d)
	{
	  while ((dir = readdir(d)) != NULL)
	  {
	    string sFilename( dir->d_name );
	    if( matchPattern( sPattern, sFilename) )
		    vsFileList.push_back( sFilename );
	  }

	  closedir(d);
	}
	
	resetListIterator();
	iFileCount = vsFileList.size();
}

int DirectoryStructure::getFileCount()
{
	return iFileCount;
}

void DirectoryStructure::getFiles( vector<string> & vsList )
{
//	vsList.clear();
//	copy( vsFileList.begin(), vsFileList.end(), vsList.begin() );
	for( unsigned int i = 0; i < vsFileList.size(); ++i)
		vsList.push_back( vsFileList[i] );
}

void DirectoryStructure::resetListIterator()
{
	itFileListItr = vsFileList.begin();
}

string DirectoryStructure::getNextFile()
{
	if ( itFileListItr != vsFileList.end() )
	{
		string sReturn = *itFileListItr;
		itFileListItr++;
		return sReturn;
	}
	else
	{
		return "";
	}
}

bool DirectoryStructure::matchPattern( const string & sPattern, const string & sFilename )
{
	// if this filename is ".", ".." or "lost+found"
	// it is not considered
	if( ( sFilename.compare( "." ) == 0 )||
		( sFilename.compare( ".." ) == 0 ) ||
		( sFilename.compare("lost+found") == 0 ) ||
		( sFilename.length() <= sPattern.size() ) )
	{
		return false;
	}
	
	// if the pattern is empty or "*", 
	// then all filenames are matched
	if( sPattern == "" || sPattern == "*" )
		return true;

	// if the filename ends with sPattern
	// then it is matched
	if ( ( sFilename.compare( sFilename.length() - sPattern.length(), sPattern.length(), sPattern ) == 0 ) )
		return true;
	else
		return false;

}

