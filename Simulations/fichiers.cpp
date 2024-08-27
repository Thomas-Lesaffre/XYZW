// Functions to open input and output files,
// read parameter values and write them in output file

#include "header.h"
#include <iostream>
#include <fstream>
#include <sstream>

using namespace std;

extern FILE * fileP;

//Opens input files:

// This opens a stream to the file defined under the token "filePar" in the header.h file.
void openfileP()
{
	fileP = fopen(filePar,"r");
	if (!fileP)
		cout << "The file " << filePar << " doesn't exist!"  << endl;
}

// The function readpar reads parameter values from input file 'par'.
/* 
The function returns 1 if it reaches the end of input file, else returns 0. It reads the parameter file character by character until it encounters an End Of File (EOF) statement, or a new line indicator (*).
It stores each of the parameter values (that have to be written in order in the 'par' file) in the parameters structure supplied as argument.
*/

bool readpar(parameters &parr)
{
	int z;
	bool term;
	do {z = fgetc(fileP);} while (!((z == '*') || (z == EOF)));
		// Lines with parameter sets must begin with *
	if (z == EOF)
	{
		cout << "\nEnd of input file\n";
		term = true;
	}
	else
	{
	
		// Parameters
		
		// General parameters
		fscanf(fileP," %lf",&parr.n);	
		fscanf(fileP," %d",&parr.tfinal);
		fscanf(fileP," %d",&parr.trig);
		fscanf(fileP," %d",&parr.eco);	


		// Initial
		fscanf(fileP," %lf",&parr.a0);					
				
		// Female
		fscanf(fileP," %lf",&parr.f0);			
		fscanf(fileP," %lf",&parr.gamma_f0);	
		fscanf(fileP," %lf",&parr.gamma_f1);	
				
		// Male	
		fscanf(fileP," %lf",&parr.m0);			
		fscanf(fileP," %lf",&parr.gamma_m0);											
		fscanf(fileP," %lf",&parr.gamma_m1);											
				
		// Genetic parameters
		fscanf(fileP," %d",&parr.nloc);					
		fscanf(fileP," %lf",&parr.u);	
		fscanf(fileP," %lf",&parr.uprom);		
		fscanf(fileP," %lf",&parr.sig);	
		
		// Repeats and measurement timings
		fscanf(fileP," %d",&parr.tmes);	 
		fscanf(fileP," %d",&parr.nmes);	 		
		fscanf(fileP," %d",&parr.n_it);	
								
																					
        term = false;
	}
	
	return term;
}

