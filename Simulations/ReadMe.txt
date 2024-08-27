This file contains the instructions to compile and run the C++ program used to produce the simulation results presented in the article "An explanation for the prevalence of XY over ZW sex determination in species derived from hermaphroditism" published in PNAS by Thomas Lesaffre, John R. Pannell and Charles Mullon.

The simulation program provided here corresponds to the version described in Appendix E, as it is the most complete, but can be used to generate all our results:
- To produce single locus simulations with fixed dominance (Appendix B.3), set parameter 'nloc' to 1 and 'uprom' to 0 so that a single locus is modelled and promoters do not mutate.
- To produce single locus simulations with evolving dominance (Appendix C.1), set parameter 'nloc' to 1 and let uprom be a positive real number between 0 and 1.
- To produce multilocus simulations (Appendix E), set parameter 'nloc' to an integer greate than one.

Program written by Thomas Lesaffre, v. 26/02/2024 

(1) Description of the files:

	- main.cpp: main file for the program, from which parameters are read and recursions are launched.
	- recursion.cpp: file where the 'recursion' function, which performs the simulation, is written.
	- functions.cpp: file where all the functions used within the 'recursion' are constructed.
	- fichiers.cpp: file containing the functions that all reading parameters from an input file.
	- ranbin.cpp: file containing the C++ recipes for sampling from various distributions (used only for the Gaussian here).
	- mt.h: MersenneTwister random number generator.

	All relevant files are extensively annotated and commented. If the reader has any question, they are welcome to e-mail T. Lesaffre at "thomas.lesaffre@unil.ch" or "thomaslesaffre.evolbiol@gmail.com"

(2) How to compile the simulation program

	To compile this program in C++11 from a Linux distribution (Ubuntu 20.04), as we did for the work presented in the article, use the command: 'g++ -std=c++11 -o ProgramName *.cpp'

(3) How to run a simulation
	
	To run a simulation, you must first check your parameter file, which should be named 'par'. Each line is a parameter set, and it must start with an *, as shown in the provided example.
	The parameter file should be in the folder in which you wish the output the be created. Then, placing yourself in this folder, you can run the simulation from terminal by typing './relative/path/to/ProgramName'.
