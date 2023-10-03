//Author: Santiago A. Flores Roman

#include <cstring> // strlen
#include <iostream> // cout
#include <sstream> // ostringstream

#include "MC.h"
#include "tools.h"

using namespace std;

MC::MC(void){
	int i, j, k;
	ostringstream boxName;

	stats.acceptanceVol = stats.rejectionVol = stats.nVolChanges = 0;
	stats.acceptSwap = stats.rejectSwap = stats.nSwaps = 0;
	for (i=0; i<MAXBOX; i++){
		stats.widomInsertions[i] = 0;
		stats.acceptance[i] = stats.rejection[i] = stats.nDisplacements[i] = 0;
		for (j=0; j<MAXSPECIES; j++) stats.widom[i][j] = 0.;
	}
	for (i=0; i<=NBINS; i++) stats.rdf[i] = 0.;
	sim.projName = "";
	sim.rdf[0] = sim.rdf[1] = -1;
	for (i=0; i<MAXBOX; i++) sim.dr[i] = 1;
	sim.dv = 1e-3;
	sim.nEquilSets = sim.nSets = sim.nStepsPerSet = 0.;
	sim.nDispAttempts = sim.nVolAttempts = sim.nSwapAttempts = 0;
	thermoSys.temp = thermoSys.volume = thermoSys.press = -1.;
	thermoSys.nBoxes = 0;
	thermoSys.nSpecies = 0;
	thermoSys.nParts = 0;
	for (i=0; i<MAXBOX; i++){
		for (j=0; j<3; j++){
			box[i].width[j] = 0;
			box[i].PBC[j] = true;
		}
		for (j=0; j<MAXSPECIES; j++){
			box[i].fluid[j].nParts = 0;
			box[i].fluid[j].vdwPot[0] = "";
			box[i].fluid[j].epsilon[0] = box[i].fluid[j].sigma[0] = 0.;
			box[i].fluid[j].mu = 0;
			for (k=0; k<MAXPART; k++){
				box[i].fluid[j].particle[k].x = 0.;
				box[i].fluid[j].particle[k].y = 0.;
				box[i].fluid[j].particle[k].z = 0.;
				box[i].fluid[j].particle[k].energy = 0.;
			}
		}
		boxName << "Box" << i; box[i].name = boxName.str();
		box[i].geometry = "bulk";
		box[i].nLayersPerWall = 0;
		box[i].deltaLayers = 0.;
		box[i].solidDens = 0.;
		box[i].nParts = 0.;
		box[i].vdwPot = 0.;
		box[i].manyBodyE = 0.;
		box[i].pairPotE = 0.;
		box[i].boxE = 0.;
		box[i].oldEnergy = 0.;
		box[i].energy = 0.;
		box[i].maxRcut = 0.;
		box[i].fix = false;
	}
	for (i=0; i<MAXSPECIES; i++){
		fluid[i].name = "";
		fluid[i].molarMass = 0.;
		fluid[i].mu = 0.;
		for (j=0; j<MAXSPECIES; j++){
			fluid[i].vdwPot[j] = "";
			fluid[i].epsilon[j] = fluid[i].sigma[j] = fluid[i].rcut[j] = 0.;
		}
	}
}
void MC::ResetStats(void){
	int i;

	stats.acceptanceVol = stats.rejectionVol = 0;
	stats.nVolChanges = 1;
	stats.acceptSwap = stats.rejectSwap = 0;
	stats.nSwaps = 1;
	for (i=0; i<thermoSys.nSpecies; i++) stats.widomInsertions[i] = 0;
	for (i=0; i<thermoSys.nBoxes; i++){
		stats.acceptance[i] = stats.rejection[i] = 0;
		stats.nDisplacements[i] = 1;
		for (int j=0; j<thermoSys.nSpecies; j++) stats.widom[i][j] = 0.;
	}
}

