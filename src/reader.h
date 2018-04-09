#pragma once

#include "entities.h"
#include <string>
#include <iostream>
#include <fstream>

int countEntriesOnFirstLine(std::string filename, const char delimiter)
{
	std::ifstream myfile(filename.c_str());
	int entries=0;
	std::string myline;
	if( getline( myfile, myline ) )
	{
		std::stringstream stream( myline );
		std::string dummy;
		while( getline(stream, dummy, delimiter) ) ++entries;
	}
	myfile.close();
	return entries;
}

int countLines(std::string filename)
{
	std::ifstream myfile(filename.c_str());
	int lines = 0;
	if (myfile.is_open())
	{
		std::string myline;
		while( getline(myfile,myline) )
		{
			++lines;
		}
	}
	myfile.close();
	return lines;
}

void fillArrayWithEntries(dtype* array, int limit, std::string line, const char delimiter)
{
	std::stringstream ss(line);
	std::string entry;
	int entries = 0;
	while( getline(ss, entry, delimiter) )
	{
		if(entries > limit) printf("\tError! entries exceeded specified limit!\n");
		array[entries] = atof(entry.c_str());
		++entries;
	}
}

DenseMatrix* readMatrix(std::string filename,const char delimiter)
{//Returns a DenseMatrix with the colmajor filled.
	printf("Reading: %s\n",filename.c_str());

	int totalLines = countLines(filename);
	printf("\tTotal number of lines in %s: %d\n",filename.c_str(),totalLines);
	
	int entriesOnFirstLine = countEntriesOnFirstLine(filename,delimiter);
	printf("\tEntries per line: %d\n",entriesOnFirstLine);

	DenseMatrix* A;

	std::ifstream myfile(filename.c_str());
	if ( myfile.is_open() and totalLines > 0 and entriesOnFirstLine > 0 )
	{
		dtype** colmajor = new dtype*[totalLines];
		for(int col=0;col<totalLines;++col) colmajor[col] = new dtype[entriesOnFirstLine];
		std::string myline;
		int linecount = 0;
		while( getline(myfile,myline) )
		{
			fillArrayWithEntries(colmajor[linecount],entriesOnFirstLine,myline,delimiter);
			++linecount;
			// std::cout << "Read line " << linecount << std::endl;
		}
		int rows = entriesOnFirstLine;
		int cols = totalLines;
		A = new DenseMatrix(colmajor,rows,cols);
	}
	else
	{
		A = NULL;
	}
	myfile.close();
	printf("Fnished reading file.\n");
	return A;
}
