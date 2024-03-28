//Author: Santiago A. Flores Roman

#include <cstring> // strlen
#include <iostream> // cout
#include <sstream> // ostringstream
#include <filesystem> // filesystem
#include <fstream> // ofstream
#include <cmath> // pi
#include "Eigen/Dense"

#include "MC.h"
#include "tools.h"

using namespace std;

void MC::PrintParams(void){
	int i, j, ithSpecies, jthSpecies;
	double maxRcut = 0.;

	cout << "Simulation parameters:" << endl;
	cout << "\tEquilibration sets: " << sim.nEquilSets << endl;
	cout << "\tProduction sets: " << sim.nSets << endl;
	cout << "\tSteps per set: " << sim.nStepsPerSet << endl;
	cout << "\tPrint every n sets: " << sim.printEvery << endl;
	cout << "\tNum. of displacement attempts: " << sim.nDispAttempts << endl;
	cout << "\tNum. of change volume attempts: " << sim.nVolAttempts << endl;
	cout << "\tNum. of swap attempts: " << sim.nSwapAttempts << endl;
	cout << endl;

	cout << "System parameters:" << endl;
	cout << "\tNum. of particles: " << thermoSys.nParts << endl;
	cout << "\tVolume (AA^3): " << thermoSys.volume << endl;
	cout << "\tExternal temperature (K): " << thermoSys.temp << endl;
	if (thermoSys.press >= 0) cout << "\tExternal pressure (Pa): " << thermoSys.press << endl;
	cout << endl;

	cout << "Fluid parameters:" << endl;
	for (i=0; i<thermoSys.nSpecies; i++){
		cout << "Species " << i << ": " << fluid(i).name << endl;
		cout << "\tMolar mass (g/mol): " << fluid(i).molarMass << endl;
		if (fluid(i).mu != 0) cout << "\tChemical potential (K): " << fluid(i).mu << endl;
	}
	cout << endl;

	cout << "Fluid-fluid interaction parameters:" << endl;
	for (i=0; i<thermoSys.nSpecies; i++){
		for (j=i; j<thermoSys.nSpecies; j++){
			cout << fluid(i).name << "-" << fluid(j).name << endl;
			cout << "\tVdW potential: " << fluid(i).vdwPot(j) << endl;
			cout << "\tsigma (AA): " << fluid(i).sigma(j) << endl;
			cout << "\tepsilon (K): " << fluid(i).epsilon(j) << endl;
			cout << "\trcut (AA): " << fluid(i).rcut(j) << endl;
		}
	}
	cout << endl;

	cout << "Box parameters:" << endl;
	for (i=0; i<thermoSys.nBoxes; i++){
		cout << "Box " << i << ": " << box(i).name << endl;
		cout << "\tGeometry: " << box(i).geometry << endl;
		cout << "\tInitial width (AA): " << box(i).width(2) << endl;
		if (box(i).solidDens > 0) cout << "\tSurface density (AA^-2): " << box(i).solidDens << endl;
		if (box(i).nLayersPerWall > 1) cout << "\tLayers per wall: " << box(i).nLayersPerWall << endl;
		if (box(i).deltaLayers > 0) cout << "\tDistance between layers (AA): " << box(i).deltaLayers << endl;
		cout << "\tInitial number of particles: " << box(i).nParts << endl;
		if (box(i).fix) cout << "\tVolume fixed." << endl;
	}
	cout << endl;

	cout << "Box-fluid interaction parameters:" << endl;
	for (i=0; i<thermoSys.nBoxes; i++){
		for (j=0; j<thermoSys.nSpecies; j++){
			cout << box(i).name << "-" << fluid(j).name << endl;
			if (box(i).fluid(j).vdwPot(0) != "") cout << "\tVdW potential: " << box(i).fluid(j).vdwPot(0) << endl;
			if (box(i).fluid(j).sigma(0) > 0) cout << "\tsigma (AA): " << box(i).fluid(j).sigma(0) << endl;
			if (box(i).fluid(j).epsilon(0) > 0) cout << "\tepsilon (K): " << box(i).fluid(j).epsilon(0) << endl;
			cout << "\tInitial number of particles: " << box(i).fluid(j).nParts << endl;
		}
	}
	cout << endl;

	if (sim.rdf(0) > -1 && sim.rdf(1) > -1){
		ithSpecies = sim.rdf(0);
		jthSpecies = sim.rdf(1);
		cout << "Sample RDF: " << fluid(ithSpecies).name << "-" << fluid(jthSpecies).name << endl;
	}
	if (sim.printTrajectory) cout << "Print Trajectory: Yes" << endl;
	else cout << "Print Trajectory: No" << endl;
	cout << endl;

	for (i=0; i<thermoSys.nSpecies; i++){
		for (j=i; j<thermoSys.nSpecies; j++){
			if (maxRcut < fluid(i).rcut(j)) maxRcut = fluid(i).rcut(j);
		}
	}
	for (i=0; i<thermoSys.nBoxes; i++){
		if ((box(i).width(2) < 2*maxRcut) && (box(i).geometry == "bulk")){
			cout << "Error: Box size must be larger than twice the cut-off radius." << endl;
			cout << "Box: " << box(i).name << endl;
			cout << endl;
			exit(EXIT_FAILURE);
		}
	}
}
void MC::OutputDirectory(void){
	Tools tls;
	ostringstream outDirName, boxDirName, command;

	outDirName << "./" << sim.projName;
	if (! std::filesystem::is_directory(outDirName.str())){
		command << "mkdir " << outDirName.str();
		system(command.str().c_str());
		command.str(""); command.clear();
	}
	for (int i=0; i<thermoSys.nBoxes; i++){
		boxDirName << outDirName.str() << "/" << box(i).name;
		if (! filesystem::is_directory(boxDirName.str())){
			command << "mkdir " << boxDirName.str();
			system(command.str().c_str());
			command.str(""); command.clear();
		}
		boxDirName.str(""); boxDirName.clear();
	}
}
// Prints the radial distribution function into a file.
// Note: RDF is printed after simulation finishes.
void MC::PrintRDF(void){
	double constant, deltaR, lowR, highR, rho, dn_ideal;
	Eigen::Matrix<double,Eigen::Dynamic,1> gofr = Eigen::Matrix<double,Eigen::Dynamic,1>(NBINS+1);
	int ithSpecies, jthSpecies, nParts;
	Tools tls;
	ofstream rdfFile;
	ostringstream outDirName, outFileName;

	if (sim.rdf(0) > -1 && sim.rdf(1) > -1){
		ithSpecies = sim.rdf(0);
		jthSpecies = sim.rdf(1);
		deltaR = fluid(ithSpecies).rcut(jthSpecies)/NBINS;
		if (ithSpecies == jthSpecies) nParts = box(0).fluid(ithSpecies).nParts;
		else nParts = box(0).fluid(ithSpecies).nParts + box(0).fluid(jthSpecies).nParts;
		rho = nParts/box(0).volume;
		constant = (4*pi*rho/3);
		for (int bin=1; bin<=NBINS; bin++){
			gofr(bin) = stats.rdf(bin)/(1.*nParts*sim.nStepsPerSet*sim.nSets);
			lowR = bin*deltaR;
			highR = lowR + deltaR;
			dn_ideal = constant*(tls.Pow(highR,3)-tls.Pow(lowR,3));
			gofr(bin) /= dn_ideal;
		}
		//Write into the file.
		outDirName << "./" << sim.projName;
		outFileName << outDirName.str() << "/" << box(0).name << "/rdf_";
		outFileName << fluid(ithSpecies).name << "-" << fluid(jthSpecies).name << ".dat";
		rdfFile.open(outFileName.str());
		rdfFile << "nBins: "<< NBINS << endl;
		rdfFile << "r(AA)\tg(r)" << endl;
		for (int bin=1; bin<=NBINS; bin++) rdfFile << (bin+0.5)*deltaR << "\t" << gofr(bin) << endl;
		rdfFile.close();
		outFileName.str(""); outFileName.clear();
		rdfFile.clear();
	}
}
void MC::PrintTrajectory(size_t set){
	Tools tls;
	ofstream trajFile;
	Particle part;
	ostringstream outDirName, outFileName;
	double tmp=0.0;

	outDirName << "./" << sim.projName;
	for (int i=0; i<thermoSys.nBoxes; i++){
		outFileName << outDirName.str() << "/" << box(i).name << "/trajectory.exyz";
		if (set == 0 || !sim.printTrajectory) trajFile.open(outFileName.str());
		else trajFile.open(outFileName.str(), ios::app);
		trajFile << box(i).nParts << "\n";
		trajFile << "Lattice=\" " << box(i).width(0) << " " << tmp << " " << tmp << " ";
		trajFile << tmp << " " << box(i).width(1) << " " << tmp << " ";
		trajFile << tmp << " " << tmp << " " << box(i).width(2) << "\" ";
		trajFile << "Properties=species:S:1:pos:R:3 "
				<< "Set: " << set << " Disp: " << sim.dr(i) <<  " VolDisp: " << sim.dv(i) << "\n";
		for (int j=0; j<thermoSys.nSpecies; j++){
			for (int k=1; k<=box(i).fluid(j).nParts; k++){
				part = box(i).fluid(j).particle(k);
				trajFile << fluid(j).name << "\t";
				trajFile << part.x << "\t" << part.y << "\t" << part.z << "\n"; //AA
			}
		}
		trajFile.close();
		outFileName.str(""); outFileName.clear();
		trajFile.clear();
	}
}
void MC::PrintLog(size_t set){
	Tools tls;
	ofstream logFile;
	ostringstream outDirName, outBoxName, outFileName;
	double density, volume, ffEnergy, energy;

	outDirName << "./" << sim.projName;
	for (int i=0; i<thermoSys.nBoxes; i++){
		for (int j=0; j<thermoSys.nSpecies; j++){
			outFileName << outDirName.str() << "/" << box(i).name;
			outFileName << "/simulation_" << fluid(j).name << ".log";
			volume = box(i).volume; //AA^3
			ffEnergy = box(i).manyBodyE + 0.5*box(i).pairPotE;
			energy = box(i).manyBodyE + 0.5*box(i).pairPotE + box(i).boxE;
			if (set == 0){
				logFile.open(outFileName.str());
				logFile << "Set\tTemp(K)\twidth(AA)\tVolume(AA^3)\t";
				logFile << "ffE/particle(K)\tsfE/particle(K)\tE/particle(K)\t";
				logFile << "NParts\t" << "Dens(g/cm^3)\t" << "muEx(K)\t" << "mu(K)\n";
			}else{
				density = box(i).fluid(j).nParts/volume; // AA^-3
				density *= fluid(j).molarMass/na*1e24; // g/cm^3
				logFile.open(outFileName.str(), ios::app);
				logFile << fixed;
				logFile << setprecision(7);
				logFile << set << "\t";
				logFile << thermoSys.temp << "\t" << box(i).width(2) << "\t" << volume << "\t";
				logFile << ffEnergy/box(i).nParts << "\t";
				logFile << box(i).boxE/box(i).nParts << "\t";
				logFile << energy/box(i).nParts << "\t";
				logFile << box(i).nParts << "\t";
				logFile << density << "\t"; // g/cm^3
				logFile	<< box(i).fluid(j).muEx << "\t"; // K
				logFile	<< box(i).fluid(j).mu << "\n"; // K
			}
			logFile.close();
			outFileName.str(""); outFileName.clear();
			logFile.clear();
		}
	}
}
void MC::PrintStats(size_t set){
	cout << fixed;
	cout << setprecision(7);
	if (set < sim.nEquilSets) cout << "Equilibrium set: " << set << endl;
	else cout << "Set: " << set << endl;
	for (int i=0; i<thermoSys.nBoxes; i++){
		cout << "Box " << box(i).name << ":" << endl;
		if (sim.nVolAttempts > 0){
			cout << "\tAcceptVolRatio; RejectVolRatio:\t";
			cout << stats.acceptanceVol(i)*1./stats.nVolChanges << "; ";
			cout << stats.rejectionVol(i)*1./stats.nVolChanges << endl;
			if (set < sim.nEquilSets) cout << "\tVolume step size: " << sim.dv(i) << endl;
		}
		if (sim.nSwapAttempts > 0){
			cout << "\tAcceptSwapRatio; RejectSwapRatio:\t";
			cout << stats.acceptSwap*1./stats.nSwaps << "; ";
			cout << stats.rejectSwap*1./stats.nSwaps << endl;
		}
		if (stats.nDisplacements(i) > 0){
			cout << "\tAcceptDispRatio; RegectDispRatio:\t";
			cout << stats.acceptance(i)*1./stats.nDisplacements(i) << "; ";
			cout << stats.rejection(i)*1./stats.nDisplacements(i) << endl;
			if (set < sim.nEquilSets) cout << "\tDisplacement step size: " << sim.dr(i) << endl;
		}
		cout << "\tNumParticles; BoxSize; Energy/Particle:\t";
		cout << box(i).nParts << "; ";
		cout << box(i).width(2) << "; ";
		cout << box(i).energy/box(i).nParts << endl;
	}
	cout << endl;
}

