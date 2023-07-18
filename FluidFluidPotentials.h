#ifndef _FLUIDFLUIDPOTENTIALS_H_
#define _FLUIDFLUIDPOTENTIALS_H_

#include <cmath>

double LJ_Energy(int i){
	double rij;
	double eps=fluid.epsilon, sig=fluid.sigma;
	double energy=0;

	for (int j=1; j<=fluid.nParts; j++){
		if (i != j){
			rij = Distance(i,j);
			if (rij < fluid.rcut) energy += 4*eps*(pow(sig/rij,12)-pow(sig/rij,6));
		}
	}
	return 0.5*energy;
}

#endif /* _FLUIDFLUIDPOTENTIALS_H_ */
