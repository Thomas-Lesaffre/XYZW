#include "header.h"
#include "mt.h"
#include <fstream>
#include <iostream>
#include <cmath>
#include <cstdlib>
#include <sstream>
#include <string>
#include <sys/stat.h>
#include <future>
#include <thread>
#include <algorithm>
#include <iterator>

using namespace std;

MTRand eng; // This is the random number generator used in the program.

// Pointers on input and output files:

FILE * fileP;

int main()
{ 	
    cout << "Program initialization" << "\n";

	parameters par; // Define the structure that will contain parameter values
	
	openfileP(); // Open the stream to the parameter file
	
	bool end { false };
	int k;
	
	end = readpar(par); // Read parameter values
	
	// Repeat the recursion par.n_it times
	for(k=0; k<par.n_it; k++)
	{
    	recursion(par, k+1);
	}	
		
	return 0;
}


