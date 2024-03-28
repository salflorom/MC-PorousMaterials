//Author: Santiago A. Flores Roman

#include <cmath> // cbrt, sqrt
#include <string> // string
#include <cstring> // strlen
#include <iostream> // cout, stod
#include <fstream> // ifstream
#include <cstdlib> // stod, stoi, strtol
#include <regex> // regex
#include <vector>
#include "Eigen/Dense"

#include "MC.h"
#include "tools.h"

using namespace std;

// Computes the length of the box along the z axis.
double MC::ComputeBoxWidth(Box& ithBox, double volume){
	cout << ""; // To avoid Core dump.
	if (ithBox.geometry == "sphere") return 2*cbrt(3*volume/(4*pi));
	else if (ithBox.geometry == "cylinder") return 2*sqrt(volume/(pi*ithBox.width(0)));
	else if (ithBox.geometry == "slit") return volume/(ithBox.width(0)*ithBox.width(1));
	else return cbrt(volume); //Bulk phase.
}
double MC::ComputeVolume(Box& ithBox){
	Tools tls;
	string geometry = ithBox.geometry;
	double xLength=ithBox.width(0), yLength=ithBox.width(1), zLength=ithBox.width(2);

	if (geometry == "sphere") return (pi/6.)*tls.Pow(zLength,3);
	else if (geometry == "cylinder") return (pi/4.)*tls.Pow(zLength,2)*xLength;
	else if (geometry == "slit") return xLength*yLength*zLength;
	else return tls.Pow(zLength,3); // Bulk phase.
}
size_t MC::ReadInputFile(string inFileName){
	int i, j, ithBox, ithSpecies, jthSpecies;
	double maxRcut;
	Eigen::Matrix<string,50,1> commands;
	string line;
	Tools tls;
	ifstream inFile;
	size_t set = 0;

	inFile.open(inFileName);
	if (!inFile.is_open()){ //Check if file was opened successfully.
		cout << "Error opening file: " << inFileName << endl;
		exit(EXIT_FAILURE);
	}
	while (getline(inFile, line)) {
		commands = tls.SplitString(line, ' ');
		for (i=0; i<4; i++) commands(i) = tls.LowerCase(commands(i));
		// Simulation parameters.
		if (commands(0) == "projectname") sim.projName = commands(1);
		else if (commands(0) == "continueaftercrash") sim.continueAfterCrash = true;
		else if (commands(0) == "productionsets") sim.nSets = stod(commands(1));
		else if (commands(0) == "equilibriumsets") sim.nEquilSets = stod(commands(1));
		else if (commands(0) == "stepsperset") sim.nStepsPerSet = stod(commands(1));
		else if (commands(0) == "printeverynsets") sim.printEvery = stod(commands(1));
		else if (commands(0) == "ndisplacementattempts") sim.nDispAttempts = stoi(commands(1));
		else if (commands(0) == "nvolumeattempts") sim.nVolAttempts = stoi(commands(1));
		else if (commands(0) == "nswapattempts") sim.nSwapAttempts = stoi(commands(1));
		else if (commands(0) == "externaltemperature") thermoSys.temp = stod(commands(1)); // K
		else if (commands(0) == "externalpressure") thermoSys.press = stod(commands(1)); // Pa
		// Species parameters.
		else if (commands(0) == "fluidname"){
			thermoSys.nSpecies++;
			fluid(thermoSys.nSpecies-1).name = commands(1);
		}else if (commands(0) == "molarmass") {fluid(thermoSys.nSpecies-1).molarMass = stod(commands(1)); // g/mol
		}else if (commands(0) == "chemicalpotential") {fluid(thermoSys.nSpecies-1).mu = stod(commands(1));} // K
		// Box parameters.
		else if (commands(0) == "boxname"){
			thermoSys.nBoxes++;
			box(thermoSys.nBoxes-1).name = commands(1);
		}else if (commands(0) == "geometry") {box(thermoSys.nBoxes-1).geometry = commands(1);
		}else if (commands(0) == "surfacedensity") {box(thermoSys.nBoxes-1).solidDens = stod(commands(1)); // AA^-2
		}else if (commands(0) == "nlayersperwall") {box(thermoSys.nBoxes-1).nLayersPerWall = stoi(commands(1));
		}else if (commands(0) == "distancebetweenlayers") {box(thermoSys.nBoxes-1).deltaLayers = stod(commands(1));
		}else if (commands(0) == "size") {box(thermoSys.nBoxes-1).width(2) = stod(commands(1)); // AA. width(2) (z axis) always keeps for pore size.
		}else if (commands(0) == "fixvolume") {box(thermoSys.nBoxes-1).fix = true;} // AA. Used only for Gibbs ensemble.
		// Interaction parameters.
		else if (commands(2) == "vdwpotential"){
			ithBox = tls.FindIndex(box, thermoSys.nBoxes, commands(0));
			ithSpecies = tls.FindIndex(fluid, thermoSys.nSpecies, commands(0));
			jthSpecies = tls.FindIndex(fluid, thermoSys.nSpecies, commands(1));
			if (ithSpecies >= 0 && jthSpecies >= 0) fluid(ithSpecies).vdwPot(jthSpecies) = commands(3);
			if (ithSpecies >= 0 && jthSpecies >= 0) fluid(jthSpecies).vdwPot(ithSpecies) = commands(3);
			if (ithBox >= 0 && jthSpecies >= 0) box(ithBox).fluid(jthSpecies).vdwPot(0) = commands(3);
		}else if (commands(2) == "sigma"){
			ithBox = tls.FindIndex(box, thermoSys.nBoxes, commands(0));
			ithSpecies = tls.FindIndex(fluid, thermoSys.nSpecies, commands(0));
			jthSpecies = tls.FindIndex(fluid, thermoSys.nSpecies, commands(1));
			if (ithSpecies >= 0 && jthSpecies >= 0) fluid(ithSpecies).sigma(jthSpecies) = stod(commands(3)); // AA
			if (ithSpecies >= 0 && jthSpecies >= 0) fluid(jthSpecies).sigma(ithSpecies) = stod(commands(3)); // AA
			if (ithBox >= 0 && jthSpecies >= 0) box(ithBox).fluid(jthSpecies).sigma(0) = stod(commands(3)); // AA
		}else if (commands(2) == "epsilon"){
			ithBox = tls.FindIndex(box, thermoSys.nBoxes, commands(0));
			ithSpecies = tls.FindIndex(fluid, thermoSys.nSpecies, commands(0));
			jthSpecies = tls.FindIndex(fluid, thermoSys.nSpecies, commands(1));
			if (ithSpecies >= 0 && jthSpecies >= 0) fluid(ithSpecies).epsilon(jthSpecies) = stod(commands(3)); // K
			if (ithSpecies >= 0 && jthSpecies >= 0) fluid(jthSpecies).epsilon(ithSpecies) = stod(commands(3)); // K
			if (ithBox >= 0 && jthSpecies >= 0) box(ithBox).fluid(jthSpecies).epsilon(0) = stod(commands(3)); // AA
		}else if (commands(2) == "rcut"){
			ithSpecies = tls.FindIndex(fluid, thermoSys.nSpecies, commands(0));
			jthSpecies = tls.FindIndex(fluid, thermoSys.nSpecies, commands(1));
			if (ithSpecies >= 0 && jthSpecies >= 0) fluid(ithSpecies).rcut(jthSpecies) = stod(commands(3)); // AA
			if (ithSpecies >= 0 && jthSpecies >= 0) fluid(jthSpecies).rcut(ithSpecies) = stod(commands(3)); // AA
		}else if (commands(2) == "numberofparticles"){
			ithBox = tls.FindIndex(box, thermoSys.nBoxes, commands(0));
			jthSpecies = tls.FindIndex(fluid, thermoSys.nSpecies, commands(1));
			if (ithBox >= 0 && jthSpecies >= 0) box(ithBox).fluid(jthSpecies).nParts = stoi(commands(3));
		}
		// Sampling parameters.
		else if (commands(0) == "samplerdf"){
			ithSpecies = tls.FindIndex(fluid, thermoSys.nSpecies, commands(1));
			jthSpecies = tls.FindIndex(fluid, thermoSys.nSpecies, commands(2));
			sim.rdf(0) = ithSpecies;
			sim.rdf(1) = jthSpecies;
		}else if (commands(0) == "printtrajectory") sim.printTrajectory = true;
		for (i=0; i<4; i++) commands(i) = "";
	}
	inFile.close();
	// Get remaining box widths and initial volumes according to the geometries given.
	maxRcut = 0.;
	for (i=0; i<thermoSys.nSpecies; i++){
		for (j=i; j<thermoSys.nSpecies; j++) if (maxRcut < fluid(i).rcut(j)) maxRcut = fluid(i).rcut(j);
	}
	for (i=0; i<thermoSys.nBoxes; i++){
		box(i).maxRcut = maxRcut;
		if (box(i).geometry == "sphere") {box(i).width(0) = box(i).width(1) = box(i).width(2);
		}else if (box(i).geometry == "cylinder"){
			box(i).width(0) = 2*maxRcut; // PBC along x axis.
			box(i).width(1) = box(i).width(2);
		}else if (box(i).geometry == "slit") {box(i).width(0) = box(i).width(1) = 2*maxRcut; // PBC along xy plane.
		}else box(i).width(0) = box(i).width(1) = box(i).width(2); // Bulk phase.
		box(i).volume = ComputeVolume(box(i));
		thermoSys.volume += box(i).volume;
	}
	// Get total number of particles and density per box.
	for (i=0; i<thermoSys.nBoxes; i++){
		for (j=0; j<thermoSys.nSpecies; j++) box(i).nParts += box(i).fluid(j).nParts;
		thermoSys.nParts += box(i).nParts;
	}
	// Set PBC restrictions.
	for (i=0; i<thermoSys.nBoxes; i++){
		if (box(i).geometry == "sphere"){
			box(i).PBC(0) = false;
			box(i).PBC(1) = false;
			box(i).PBC(2) = false;
		}else if (box(i).geometry == "cylinder"){
			box(i).PBC(0) = true;
			box(i).PBC(1) = false;
			box(i).PBC(2) = false;
		}else if (box(i).geometry == "slit"){
			box(i).PBC(0) = true;
			box(i).PBC(1) = true;
			box(i).PBC(2) = false;
		}else if (box(i).geometry == "bulk"){
			box(i).PBC(0) = true;
			box(i).PBC(1) = true;
			box(i).PBC(2) = true;
		} }
	// Set chemical potential to boxes if given.
	for (i=0; i<thermoSys.nBoxes; i++){
		for (j=0; j< thermoSys.nSpecies; j++) box(i).fluid(j).mu = fluid(j).mu;
	}
	// Set step size.
	for (i=0; i<thermoSys.nBoxes; i++) sim.dr(i) = 0.1*box(i).width(2); //AA
	// Set parameters from last configuration saved if user requested to continue after crash.
	if (sim.continueAfterCrash) set = ReadLogFile();
	return set;
}
size_t MC::ReadLogFile(void){
	Tools tls;
	string line;
	Eigen::Matrix<string,Eigen::Dynamic,1> cells;
	ostringstream inDirName, simFileName;
	ifstream simFile;
	int lineIdx;
	size_t set=0;

	inDirName << "./" << sim.projName;
	for (int i=0; i<thermoSys.nBoxes; i++){
		for (int j=0; j<thermoSys.nSpecies; j++){
			simFileName << inDirName.str() << "/" << box(i).name << "/simulation_" << fluid(j).name << ".log";
			cout << "Reading log file: " << simFileName.str() << endl;
			simFile.open(simFileName.str());
			if (simFile.is_open()){ //Check if file was opened successfully.
				lineIdx = 0;
				while (getline(simFile, line)) {
					if (lineIdx < 1){lineIdx++; continue;} // Avoid header.
					cells = tls.SplitString(line, '\t');
					set = stol(cells(0));
					thermoSys.temp = stod(cells(1));
					box(i).width(2) = stod(cells(2));
					box(i).boxE = stod(cells(5));
					box(i).fluid(j).nParts = stoi(cells(7));
				}
			}else{
				cout << "\tWarning: Last configuration not found." << endl;
				cout << "\t\tCurrent set, temperature, box width, and num. of particles will be assigned according to the input file." << endl;
			}
			cout << "\tFinished reading log file." << endl << endl;
			simFileName.str(string());
			simFileName.clear();
			simFile.close();
			simFile.clear();
		}
		box(i).nParts = 0;
		for (int j=0; j<thermoSys.nSpecies; j++) box(i).nParts += box(i).fluid(j).nParts;
	}
	thermoSys.nParts = 0;
	for (int i=0; i<thermoSys.nBoxes; i++) thermoSys.nParts += box(i).nParts;
	cout << "Continue after crash: Yes" << endl;
	cout << "Current system state:" << endl;
	cout << "\tCurrent set: " << set << endl;
	cout << "\tSystem temperature: " << thermoSys.temp << endl;
	for (int i=0; i<thermoSys.nBoxes; i++){
		cout << "\t" << box(i).name << " size: " << box(i).width(2) << endl;
		for (int j=0; j<thermoSys.nSpecies ; j++){
			cout << "\t\tNum. of particles of species " << fluid(j).name << " in the box: " << box(i).fluid(j).nParts << endl;
		}
		cout << "\tNum. of particles in the box: " << box(i).nParts << endl;
	}
	cout << endl;
	return set;
}
void MC::ReadTrajectory(void){
	Tools tls;
	string line;
	Eigen::Matrix<string,Eigen::Dynamic,1> cells;
	ostringstream inDirName, simFileName;
	ifstream simFile;
	int lineIdx=0, ithPart=0, nParts=0;

	inDirName << "./" << sim.projName;
	for (int i=0; i<thermoSys.nBoxes; i++){
		simFileName << inDirName.str() << "/" << box(i).name << "/trajectory.exyz";
		cout << "Reading trajectory file: " << simFileName.str() << endl;
		simFile.open(simFileName.str());
		if (simFile.is_open()){ //Check if file was opened successfully.
			while (getline(simFile, line)){
				// Check if there is a folowing configurations.
				if (lineIdx-2 == nParts) lineIdx = ithPart = 0;
				if (lineIdx == 0){
					cells = tls.SplitString(line, '\t');
					nParts = stoi(cells(0));
				}
				if (lineIdx == 1){
					cells = tls.SplitString(line, ' ');
					sim.dr(i) = stod(cells(14));
					sim.dv(i) = stod(cells(16));
				}
				// Read a configuration registered in the trajectory file.
				if (lineIdx < 2){lineIdx++; continue;}
				cells = tls.SplitString(line, '\t');
				for (int j=0; j<thermoSys.nSpecies; j++){
					if (cells(0) == fluid(j).name){
						box(i).fluid(j).particle(ithPart).x = stod(cells(1));
						box(i).fluid(j).particle(ithPart).y = stod(cells(2));
						box(i).fluid(j).particle(ithPart).z = stod(cells(3));
					}
				}
				lineIdx++;
				ithPart++;
			}
		}else{
			cout << "\tWarning: Last configuration not found." << endl;
			cout << "\t\tConfiguration will be created by assigning random positions." << endl;
			sim.continueAfterCrash = false;
			InitialConfig();
			sim.continueAfterCrash = true;
		}
		cout << "\tFinished reading configuration." << endl;
		simFileName.str(string());
		simFileName.clear();
		simFile.close();
		simFile.clear();
	}
	cout << endl;
}

