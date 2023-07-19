#ifndef _FLUIDFLUIDPOTENTIALS_H_
#define _FLUIDFLUIDPOTENTIALS_H_

#include <cmath>
#include <string>
#include <vector>
#include <unordered_map>

#include "MC.h"

using namespace std;

double Distance(vector<double> posi, vector<double> posj, double boxWidth){
	float dist, dx, dy, dz;

	dx = posj[0] - posi[0]; //AA. x-axis
	dx -= round(dx/boxWidth) * boxWidth; //AA
	dy = posj[1] - posi[1]; //AA. y-axis.
	dy -= round(dy/boxWidth) * boxWidth; //AA
	dz = posj[2] - posi[2]; //AA. z-axis.
	dz -= round(dz/boxWidth) * boxWidth; //AA
	dist = sqrt(dx*dx + dy*dy + dz*dz); //AA
	return dist; //AA
}
double LJ_Energy(int index, unordered_map<string,double> params, vector<vector<double>> pos){
	double rij, eps, sig, rcut, boxWidth;
	int nParts;
	double energy=0;

	eps = params["epsilon[K]"];
	sig = params["sigma[AA]"];
	rcut = params["rcut[AA]"];
	nParts = int(params["numOfParticles"]);
	boxWidth = params["boxWidth[AA]"];
	for (int j=1; j<=nParts; j++){
		if (index != j){
			rij = Distance(pos[index], pos[j], boxWidth);
			if (rij < rcut) energy += 4*eps*(pow(sig/rij,12)-pow(sig/rij,6));
		}
	}
	return 0.5*energy; //K
}

#endif /* _FLUIDFLUIDPOTENTIALS_H_ */
