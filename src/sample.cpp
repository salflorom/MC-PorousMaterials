//Author: Santiago A. Flores Roman

#include <cmath> // pi, exp
#include <cstring> // strlen

#include "MC.h"
#include "tools.h"

using namespace std;

// Computes mu in each box from the Widom insearsions performed in the ComputeWidom method.
// Note: It only works for simple fluids.
// Note 2: the chemical potential is printed in all boxes.
void MC::ComputeChemicalPotential(void){
	Tools tls;
	double thermalWL, particleMass, volume, muIdeal, muExcess, insertionParam;
	double extPress;

	if (thermoSys.nBoxes == 1 && sim.nSwapAttempts == 0){
		particleMass = fluid[0].molarMass/na*1e-3; // kg/particle
		thermalWL = planck/sqrt(2.0*pi*particleMass*kb*thermoSys.temp)*1e10; //AA
		volume = box[0].volume; //AA^3
		insertionParam = stats.widom[0][0]/stats.widomInsertions;
		if (sim.nVolAttempts > 0){
			extPress = thermoSys.press/kb*1e-30; // K/AA^3
			muIdeal = thermoSys.temp*log(tls.Pow(thermalWL,3)*extPress/thermoSys.temp); //K
			muExcess = -thermoSys.temp*log(insertionParam*extPress/thermoSys.temp); //K
		}else{
			muIdeal = thermoSys.temp*log(tls.Pow(thermalWL,3)*(box[0].fluid[0].nParts+1)/volume); //K
			muExcess = -thermoSys.temp*log(insertionParam); //K
		}
	}else if (thermoSys.nBoxes > 1){
		muIdeal = thermoSys.temp*log(tls.Pow(thermalWL,3)); //K
		for (int i=0; i<thermoSys.nSpecies; i++){
			insertionParam = stats.widom[0][i]/stats.widomInsertions;
			muExcess = -thermoSys.temp*log(insertionParam); //K
			box[0].fluid[i].muEx = muExcess; //K
			box[0].fluid[i].mu = muIdeal+muExcess; //K
			insertionParam = stats.widom[1][i]/stats.widomInsertions;
			muExcess = -thermoSys.temp*log(insertionParam); //K
			box[1].fluid[i].muEx = muExcess; //K
			box[1].fluid[i].mu = muIdeal+muExcess; //K
		}
	}
}
// Computes the radial distribution function every MC cycle.
// Note: Due to the expensive PC cost, it only computes the RDF in box 0 for the chosen pair of species.
void MC::ComputeRDF(void){
	int bin, ithSpecies, jthSpecies;
	double deltaR;
	Tools tls;
	Particle ithPart, jthPart;

	if (sim.rdf[0] > -1 && sim.rdf[1] > -1){
		ithSpecies = sim.rdf[0];
		jthSpecies = sim.rdf[1];
		deltaR = fluid[ithSpecies].rcut[jthSpecies]/(1.*NBINS);
		if (ithSpecies == jthSpecies){
			for (int i=1; i<=int(box[0].fluid[ithSpecies].nParts)-1; i++){
				for (int j=i+1; j<=int(box[0].fluid[ithSpecies].nParts); j++){
					ithPart = box[0].fluid[ithSpecies].particle[i];
					jthPart = box[0].fluid[ithSpecies].particle[j];
					bin = int(NeighDistance(0, ithPart, jthPart)/deltaR)+1;
					if (bin <= NBINS) stats.rdf[bin] += 2; // Takes into account both i->j and j->i.
				}
			}
		}else{
			for (int i=1; i<=int(box[0].fluid[ithSpecies].nParts); i++){
				for (int j=1; j<=int(box[0].fluid[jthSpecies].nParts); j++){
					ithPart = box[0].fluid[ithSpecies].particle[i];
					jthPart = box[0].fluid[jthSpecies].particle[j];
					bin = int(NeighDistance(0, ithPart, jthPart)/deltaR)+1;
					if (bin <= NBINS) stats.rdf[bin] += 2; // Takes into account both i->j and j->i.
				}
			}
		}
	}
}

