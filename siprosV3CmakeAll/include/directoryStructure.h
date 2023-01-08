
/*
 * Class used to manage files in a given directory.
 * Reads a particular directory and gets all the *.chro files from that
 * directory. Has utility functions to provide file count and individual
 * files.
 */

#ifndef DIRECTORYSTRUCTURE_H
#define DIRECTORYSTRUCTURE_H



#ifndef _WIN32
#include <unistd.h>
#endif


#include <cstring>
#include <iostream>
#include <vector>
#include <string>
#include <sys/types.h>
#include <dirent.h>

using namespace std;

class DirectoryStructure
{
	public:
		/*
		 * Create an instance by providing a directory name.
		 */
		DirectoryStructure();
		DirectoryStructure( const string & );

		virtual ~DirectoryStructure();

		/* 
		 * Set the directory explicitly, use it when you want to use
		 * the same instance for a different directory.
		 */
		void 	setDirectory( const string & );
		
		/*
		 * Set the tail pattern for the files to be detected.
		 * Eg: To detect *.chro files, set the pattern as
		 * ".chro".
		 *
		 * BEWARE: This class doesn't support full fledged pattern
		 * matching
		 */
		virtual void setPattern( string sPattern );

		/*
		 * Get the number of matched files 
		 */
		int  	getFileCount();

		/*
		 * Get all the files at a time.
		 */
		void 	getFiles( vector<string> & vsList );

		/*
		 * Re-entrant function, gives you serial access to all the files.
		 */
		// start from the first file 
		void	resetListIterator();
		// return the next filename,
		// return empty string "", if there is no more file
		string 	getNextFile();


	private:

		// the total number of files in this directory
		int 	iFileCount;

		// the name of the directory
		string 	sDirectory;

		// the list of file names in this directory
		vector<string> vsFileList;

		// iterator used for serial access
		vector<string>::iterator itFileListItr;


		// determine if this filename matches this pattern
		bool matchPattern(  const string & sPattern, const string & sFilename );

};

#endif //DIRECTORYSTRUCTURE_H

