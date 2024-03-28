//Author: Santiago A. Flores Roman

#include <cstring> // strlen
#include <iostream> // cout
#include <sstream> // ostringstream
#include "Eigen/Dense"

#include "MC.h"
#include "tools.h"

using namespace std;

MC::MC(void){
	int i, j, k;
	ostringstream boxName;

	stats.nSwaps = 0;
	stats.acceptSwap = stats.rejectSwap = 0;
	stats.nVolChanges = 0;
	stats.acceptanceVol = Eigen::VectorXi::Zero(MAXBOX);
	stats.rejectionVol = Eigen::VectorXi::Zero(MAXBOX);
	stats.acceptance = Eigen::VectorXi::Zero(MAXBOX);
	stats.rejection = Eigen::VectorXi::Zero(MAXBOX);
	stats.nDisplacements = Eigen::VectorXi::Zero(MAXBOX);
	stats.barC = Eigen::MatrixXd::Ones(MAXBOX,MAXSPECIES);
	stats.widomInsertions = Eigen::MatrixXi::Ones(MAXBOX,MAXSPECIES);
	stats.widomDeletions = Eigen::MatrixXi::Ones(MAXBOX,MAXSPECIES);
	stats.widomInsert = Eigen::MatrixXd::Zero(MAXBOX,MAXSPECIES);
	stats.rdf = Eigen::VectorXd::Zero(NBINS);
	sim.projName = "";
	sim.printTrajectory = sim.continueAfterCrash = false;
	sim.rdf << -1, -1;
	sim.dr = Eigen::VectorXd::Ones(MAXBOX);
	sim.dv = Eigen::VectorXd::Constant(MAXBOX,1,1e-3);
	sim.nEquilSets = sim.nSets = sim.nStepsPerSet = 0.;
	sim.nDispAttempts = sim.nVolAttempts = sim.nSwapAttempts = 0;
	thermoSys.temp = thermoSys.volume = thermoSys.press = -1.;
	thermoSys.nBoxes = 0;
	thermoSys.nSpecies = 0;
	thermoSys.nParts = 0;
	for (i=0; i<MAXBOX; i++){
		box(i).PBC = Eigen::Vector3i::Ones();
		box(i).width = Eigen::Vector3d::Zero();
		for (j=0; j<MAXSPECIES; j++){
			box(i).fluid(j).nParts = 0;
			box(i).fluid(j).vdwPot(0) = "";
			box(i).fluid(j).epsilon(0) = box(i).fluid(j).sigma(0) = 0.;
			box(i).fluid(j).mu = 0;
			for (k=0; k<MAXPART; k++){
				box(i).fluid(j).particle(k).x = 0.;
				box(i).fluid(j).particle(k).y = 0.;
				box(i).fluid(j).particle(k).z = 0.;
				box(i).fluid(j).particle(k).energy = 0.;
			}
		}
		boxName << "Box" << i; box(i).name = boxName.str();
		box(i).geometry = "bulk";
		box(i).nLayersPerWall = 0;
		box(i).deltaLayers = 0.;
		box(i).solidDens = 0.;
		box(i).nParts = 0.;
		box(i).vdwPot = 0.;
		box(i).manyBodyE = 0.;
		box(i).pairPotE = 0.;
		box(i).boxE = 0.;
		box(i).oldEnergy = 0.;
		box(i).energy = 0.;
		box(i).maxRcut = 0.;
		box(i).fix = false;
	}
	for (i=0; i<MAXSPECIES; i++){
		fluid(i).name = "";
		fluid(i).molarMass = 0.;
		fluid(i).mu = 0.;
		fluid(i).epsilon = Eigen::VectorXd::Zero(MAXSPECIES);
		fluid(i).sigma = Eigen::VectorXd::Zero(MAXSPECIES);
		fluid(i).rcut = Eigen::VectorXd::Zero(MAXSPECIES);
		for (j=0; j<MAXSPECIES; j++) fluid(i).vdwPot(j) = "";
	}
}
void MC::ResetStats(void){
	int i;

	stats.nVolChanges = 1;
	stats.acceptSwap = stats.rejectSwap = 0;
	stats.nSwaps = 1;
	for (i=0; i<thermoSys.nBoxes; i++){
		for (int j=0; j<thermoSys.nSpecies; j++){
			stats.barC(i,j) = box(i).fluid(j).muEx - thermoSys.temp*log(stats.widomInsertions(i,j)/stats.widomDeletions(i,j));
			stats.barC(i,j) = 1.;
			stats.widomInsertions(i,j) = stats.widomDeletions(i,j) = 1;
			stats.widomInsert(i,j) = stats.widomDelete(i,j) = 0.;
		}
		stats.acceptanceVol(i) = stats.rejectionVol(i) = 0;
		stats.acceptance(i) = stats.rejection(i) = 0;
		stats.nDisplacements(i) = 1;
	}
}

