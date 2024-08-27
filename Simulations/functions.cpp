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

// The function 'male' takes the genotype of an individual (two chromosomes) and parameter values, and returns its male fecundity.
double male(double xv, parameters parv, int time)
{
	double pol;
	
	if(time < parv.eco)
	{
		pol = parv.m0*pow(1 - xv, parv.gamma_m0); 	
	}else{	
		pol = parv.m0*pow(1 - xv, parv.gamma_m1); 		
	}
	
	return(pol); // Return male fecundity.
}

// Female gain curve: Exactly the same procedure as for male fecundity, but applied to the female fecundity.
double female(double xv, parameters parv, int time)
{
	double ovu;
	
	if(time < parv.eco)
	{
		ovu = parv.f0*pow(xv, parv.gamma_f0);
	}else{
		ovu = parv.f0*pow(xv, parv.gamma_f1);
	} 
	
	return(ovu); // Return female fecundity.
}


// Function used to compute individual strategies: exactly the same approach as above, but return allocation instead of fecundity.
double alloc(chrom &c1, chrom &c2, parameters parv)
{
	double x, h, c;
	int i;
	
	x = 0.0;
	c = 0.0;
	h = 0.0;
	
	for(i=0; i<parv.nloc; i++)
	{
		c = (c1.c[i] + c2.c[i])/2.0; // Contribution of the locus
		h += c; // Sum of all contributions
		x += c*( c1.m[i].exp*c1.m[i].sex + c2.m[i].exp*c2.m[i].sex )/(c1.m[i].exp + c2.m[i].exp); // Effect of the locus on the phenotype... 
	}
	
	x /= h; // ...divided by total contributions
	
	// Safety check
	if(x < 0){x = 0;}
	if(x > 1){x = 1;}

	return(x);
} 

// This function performs meiosis taking the two parental chromosomes, parameter values and time since the start of the simulation as arguments, and returning a chromosome.
chrom rec_mut(chrom &c1, chrom &c2, parameters parv, int time)
{
	// Defining some variables:
	chrom off; // Resulting offspring chromosome
	int i; // Iterator
	
	// Mendelian segregation at the contribution modifier
	if(eng.rand() < 0.5){
		off.c = c1.c;
	}else{
		off.c = c2.c;
	}
	
	// Mutation at the contribution modifier with probability 'u'
	if(eng.rand() < parv.u){
		for(i=0; i<parv.nloc; i++){
			off.c[i] += parv.sig*gaussdev();
			if(off.c[i] < 0.000001)
			{
				off.c[i] = 0.000001;
			}
		}
	}
		
	 // For each locus...
	for(i=0; i<parv.nloc; i++)
	{
		/* 
		Inherit one of the two alleles (sex allocation allele + its promoter) from the parental chromosomes with equal probability. 
		This corresponds to free recombination between the loci when parv.nloc > 1
		*/
		if(eng.rand() < 0.5) 
		{
			off.m.push_back( c1.m[i] );
		}
		else
		{
			off.m.push_back( c2.m[i] );		
		}
		
		if(eng.rand() < parv.u) // If mutation occurs at the sex allocation locus...
		{
			off.m[i].sex += parv.sig*gaussdev(); // ...mutate the sex allocation allele.			
			
			if(off.m[i].sex < 0){ off.m[i].sex = 0; }
			if(off.m[i].sex > 1){ off.m[i].sex = 1; }
		}

		if((eng.rand() < parv.uprom)&&(time > parv.trig)) // If mutation occurs at the promoter AND is allowed (time > parv.trig),
		{
			off.m[i].exp += parv.sig*gaussdev(); // Mutate promoter affinity and keep it within bounds.
			
			if(off.m[i].exp < 0.000001)
			{
				off.m[i].exp = 0.000001;
			}
		}		
	} 
	
	return(off); // Return the resulting chromosome.
}

