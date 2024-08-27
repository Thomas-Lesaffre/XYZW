#ifndef DEPRESSION_H
#define DEPRESSION_H

#include <vector>
#include <iostream>
#include <future>
#include <thread>
#include <iterator>

using namespace std;

// global variables

#define filePar "par" 

struct parameters
{
	// General parameters
	double n; // Population size
	int tfinal; // # of generations 
	int trig; // # of generations after which dominance evolution is triggered
	int eco; // # of gen before ecological change.
	
	// Initial conditions
	
	double a0;
	
	// Female
	double f0; // Female fecundity constant
	double gamma_f0;	// Exponent of the female gain curve
	double gamma_f1;
	
	// Male	
	double m0;	// Male fecundity constant		
	double gamma_m0;	// Exponent of the male gain curve										
	double gamma_m1;
	
	// Genetic parameters
	int nloc; // # of sex allocation loci (always set to 1 here)
	double u; // Mutation rate
	double uprom; // Mutation rate at the promoter.
	double sig;	// Size of mutation steps (StDev of the Gaussian)
	
	// Repeats and measurements
	int tmes; // # time between measurements.	
	int nmes; // # of measurements carried out every tmes.		
	int n_it; // # of repeats.
};

// Structure describing a "locus" (sex allocation locus and its promoter)
struct loc
{
	double sex; // Sex allocation allele carried
	double exp; // Promoter affinity of the allele carried
};

// Structure describing a chromosome as an ordered collection of loci.
struct chrom
{
	vector<loc> m; // Sex allocation loci
	vector<double> c; // Contribution modifier
};

// Prototypes of functions


void openfileP(); // Open the parameters file

bool readpar(parameters &parr); // Read parameter values

double male(double xv, parameters parv, int time); // Male gain curve
double female(double xv, parameters parv, int time); // Female gain curve
double alloc(chrom &c1, chrom &c2, parameters parv); // Function used to calculate an individual's expressed sex allocation phenotype given its genotype
chrom rec_mut(chrom &c1, chrom &c2, parameters parv, int time); // Function used to build the offspring's chromosome during reproduction.

void recursion(parameters parv, int it);

void cntl_c_handler(int bidon);

// Distributions

double gammln(const double xx);
double poisdev(const double xm);
double binldev(const double pp, const int n);
double gaussdev();

#endif
