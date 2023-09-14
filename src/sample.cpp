//Author: Santiago A. Flores Roman

#include <cmath> // log
#include <cstring> // strlen

#include "MC.h"
#include "tools.h"

using namespace std;

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

