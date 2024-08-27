#include "header.h"
#include "mt.h"
#include <vector>
#include <cmath>
#include <algorithm>
#include <fstream>
#include <sstream>
#include <future>
#include <chrono>
#include <thread>

using namespace std;

extern MTRand eng;

// Recursion function

void recursion(parameters parv, int it)
{
	
	// Creating output files for phenotypes, sex allocation genotypes and promoter affinities
	
	// Output file for phenotypes
	char nomFichier[256];
	stringstream nomF;
	nomF << "strategies_gf" << parv.gamma_f0 << "_" << parv.gamma_f1 << "_gm" << parv.gamma_m0 << "_" << parv.gamma_m1 << "_eco" << parv.eco << "-" << 	 it;	
	nomF >> nomFichier;
	ofstream fout(nomFichier); // Output stream to the file
	
	// Output file for the contribution of each locus to phenotype
	char nomFichierg[256];
	stringstream nomFg;
	nomFg << "contributions_gf" << parv.gamma_f0 << "_" << parv.gamma_f1 << "_gm" << parv.gamma_m0 << "_" << parv.gamma_m1 << "_eco" << parv.eco << "-" << it;		
	nomFg >> nomFichierg;
	ofstream foutg(nomFichierg); // Output stream to the file
	
	// Output file for the contribution of each locus to phenotype
	char nomFichiera[256];
	stringstream nomFa;
	nomFa << "affinities_gf" << parv.gamma_f0 << "_" << parv.gamma_f1 << "_gm" << parv.gamma_m0 << "_" << parv.gamma_m1 << "_eco" << parv.eco << "-" << it;		
	nomFa >> nomFichiera;
	ofstream fouta(nomFichiera); // Output stream to the file

	// Output file for the contribution of each locus to phenotype
	char nomFichiere[256];
	stringstream nomFe;
	nomFe << "effects_gf" << parv.gamma_f0 << "_" << parv.gamma_f1 << "_gm" << parv.gamma_m0 << "_" << parv.gamma_m1 << "_eco" << parv.eco << "-" << it;		
	nomFe >> nomFichiere;
	ofstream foute(nomFichiere); // Output stream to the file	
	
	// Defining some variables:
	chrom cdum; // Dummy chromosome struct
	loc ldum; // Dummy locus struct
	
	int i, j, k, t, nb, stop;
	double r, x, x0, xx0, tot_cont;

	vector<chrom> pop; // Vector containing the population
	vector<chrom> par; // Vector containing the parental population
	
	vector<double> fecm, fecf, cont(parv.nloc), ccont(parv.nloc); // Vectors of male and female fecundities

	// Initialising the population

	x0 = parv.gamma_f0/(parv.gamma_f0 + parv.gamma_m0); 
	
	/* 
	We initialise the population so that it already expresses the singular strategy here.
	This can be changed by setting another value for 'x0' and recompiling the program.
	The reason we did that was to save some time when running many replicates for the evolution of dominance.
	*/
	
	cdum.m.clear(); // Clear the dummy chromsome
	
	for(j=0; j<parv.nloc; j++) // Fill it with the appropriate number of loci
	{
		ldum.sex = x0; // Monomorphic for some sex allocation strategy...
		ldum.exp = parv.a0; // ... and affinity
					
		cdum.m.push_back(ldum); 
		cdum.c.push_back(parv.a0); // Contribution modifier is initialised for equal contributions of all loci.
	}	
	
	for(i=0; i<2*parv.n; i++) // For each chromosome (2 per individual)
	{		
		pop.push_back( cdum ); // Push the chromosome into the population and parental population vectors
		par.push_back( cdum );
	}
	
	for(t=0; t<parv.tfinal; t++) // For each time step,
	{		
		// Create a vector of cumulated male and female fecundities
		fecm.clear();
		fecf.clear();
				
		for(i=0; i<parv.n; i++) // For each individual,
		{
			nb = 2*i;
			
			// Copy it in the parental pop vector
			par[nb] = pop[nb];
			par[nb+1] = pop[nb+1];
			
			// Determined its fecundity and add it to the cumulated fecundity vector of each sex.
			if(i == 0) 
			{
				fecm.push_back( male( alloc(par[nb], par[nb+1], parv), parv, t) );			
				fecf.push_back( female( alloc(par[nb], par[nb+1], parv), parv, t) );							
			}
			else
			{
				fecm.push_back( fecm.back() + male( alloc(par[nb], par[nb+1], parv), parv, t) );			
				fecf.push_back( fecf.back() + female( alloc(par[nb], par[nb+1], parv), parv, t) );			
			}

		}
						
		// Sampling the next generation
		for(i=0; i<parv.n; i++) // For each breeding spot,
		{
			nb = 2*i;
			
			// Select a mother
			/*
			We sample a random number between 0 and the sum of all fecundities, and then go down the cumulated fecundity vector until we find the first value that is larger or equal to the sampled number.
			This way, we sample individual 'k' with a probability that depends on its fecundity (the more fecund an individual is, the more likely it will be that the sampled value falls in its range.
			*/
			k=-1;
			r=eng.rand( fecf.back() );
			do{
				k++;
			}while( fecf[k] < r );
			
			pop[nb] = rec_mut(par[2*k], par[2*k + 1], parv, t); // Create the maternally inherited chromosome
			
			// Select a father
			k=-1;
			r=eng.rand( fecm.back() );
			do{
				k++;
			}while( fecm[k] < r );
			
			pop[nb+1] = rec_mut(par[2*k], par[2*k + 1], parv, t); // Create the paternally inherited chromosome
		}
		
		if(t % parv.tmes == 0) // Every tmes generations
		{
			cout << "t=" << t << endl;		
			
			// Initialise the vector of contributions
			for(j=0; j<parv.nloc; j++){
				ccont[j] = 0.0;
			}
			
			/// For each of the nmes individuals measured
			for(i=0; i<parv.nmes; i++){
				nb=2*i;
				
				// Sex allocation expressed
				x=alloc(pop[nb],pop[nb+1],parv);
				if(i==0){
					fout << x;
				}else{
					fout << " " << x;
				}
				
				// Affinities and effects at each locus
				if(i==0){
					for(j=0; j<parv.nloc; j++){
						if(j==0){
							fouta << pop[nb].m[j].exp << " " << pop[nb+1].m[j].exp;		
							foute << pop[nb].m[j].sex << " " << pop[nb+1].m[j].sex;																									
						}else{
							fouta << " " << pop[nb].m[j].exp << " " << pop[nb+1].m[j].exp;		
							foute << " " << pop[nb].m[j].sex << " " << pop[nb+1].m[j].sex;																																							
						}
					}
				}else{
					for(j=0; j<parv.nloc; j++){
						fouta << " " << pop[nb].m[j].exp << " " << pop[nb+1].m[j].exp;
						foute << " " << pop[nb].m[j].sex << " " << pop[nb+1].m[j].sex;																																								
					}
				}
				
				// Average contributions of each locus
				tot_cont=0.0;
				for(j=0; j<parv.nloc; j++){
					tot_cont += (pop[nb].c[j] + pop[nb+1].c[j])/2.0;
					cont[j] = (pop[nb].c[j] + pop[nb+1].c[j])/2.0;
				}
				
				for(j=0; j<parv.nloc; j++){
					cont[j] /= tot_cont;
					ccont[j] += cont[j];
				}
			}
			
			for(j=0; j<parv.nloc; j++){
				ccont[j] /= parv.nmes;
				if(j==0){
					foutg << ccont[j];
				}else{
					foutg << " " << ccont[j];				
				}
			}
				
			fout << endl;
			foutg << endl;			
			fouta << endl;
			foute << endl;	
			
			// Stopping system for simulations with dominance evolution. This is allows the simulation run to dynamically stop once no individual expresses an intermediate sex allocation value (so that dioecy has evolved).
			stop=1;
			for(i=0; i<parv.n; i++){
				nb=2*i;
				if( alloc(pop[nb], pop[nb+1], parv) > 0.05 && alloc(pop[nb], pop[nb+1], parv) < 0.95  ){
					stop=0;
					break;
				}
			}
			
				
		} // End of 'if' for writing	

		if(stop == 1){
			break;
		}	
	} // End of loop over time	
} // End of function
