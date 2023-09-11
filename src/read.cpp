//Author: Santiago A. Flores Roman

#include <cmath> // pi, exp
#include <cstring> // strlen
#include <iostream> // cout

#include "MC.h"
#include "tools.h"

using namespace std;

// Computes the length of the box along the z axis.
double MC::ComputeBoxWidth(Box ithBox, double volume){
	cout << ""; // To avoid Core dump.
	if (ithBox.geometry == "sphere") return 2*cbrt(3*volume/(4*pi));
	else if (ithBox.geometry == "cylinder") return 2*sqrt(volume/(pi*ithBox.width[0]));
	else if (ithBox.geometry == "slit") return volume/(ithBox.width[0]*ithBox.width[1]);
	else return cbrt(volume); //Bulk phase.
}
double MC::ComputeVolume(Box ithBox){
	Tools tls;
	string geometry = ithBox.geometry;
	double xLength=ithBox.width[0], yLength=ithBox.width[1], zLength=ithBox.width[2];

	if (geometry == "sphere") return (pi/6.)*tls.Pow(zLength,3);
	else if (geometry == "cylinder") return (pi/4.)*tls.Pow(zLength,2)*xLength;
	else if (geometry == "slit") return xLength*yLength*zLength;
	else return tls.Pow(zLength,3); // Bulk phase.
}
void MC::ReadInputFile(string inFileName){
	int i, j, ithBox, ithSpecies, jthSpecies;
	double maxRcut;
	string* commands;
	string line;
	Tools tls;
	ifstream inFile;

	inFile.open(inFileName);
	if (!inFile.is_open()){ //Check if file was opened successfully.
		cout << "Error opening file" << endl;
		exit(EXIT_FAILURE);
	}
	while (getline(inFile, line)) {
		commands = tls.SplitString(line, ' ');
		for (i=0; i<4; i++) commands[i] = tls.LowerCase(commands[i]);
		// Simulation parameters.
		if (commands[0] == "projectname") sim.projName = commands[1];
		else if (commands[0] == "productionsets") sim.nSets = stod(commands[1]);
		else if (commands[0] == "equilibriumsets") sim.nEquilSets = stod(commands[1]);
		else if (commands[0] == "cyclesperset") sim.nCyclesPerSet = stod(commands[1]);
		else if (commands[0] == "printeverynsets") sim.printEvery = stod(commands[1]);
		else if (commands[0] == "ndisplacementattempts") sim.nDispAttempts = stoi(commands[1]);
		else if (commands[0] == "nvolumeattempts") sim.nVolAttempts = stoi(commands[1]);
		else if (commands[0] == "nswapattempts") sim.nSwapAttempts = stoi(commands[1]);
		else if (commands[0] == "externaltemperature") thermoSys.temp = stod(commands[1]); // K
		else if (commands[0] == "externalpressure") thermoSys.press = stod(commands[1]); // Pa
		// Species parameters.
		else if (commands[0] == "fluidname"){
			thermoSys.nSpecies++;
			fluid[thermoSys.nSpecies-1].name = commands[1];
		}else if (commands[0] == "molarmass") {fluid[thermoSys.nSpecies-1].molarMass = stod(commands[1]); // g/mol
		}else if (commands[0] == "chemicalpotential") {fluid[thermoSys.nSpecies-1].mu = stod(commands[1]);} // K
		// Box parameters.
		else if (commands[0] == "boxname"){
			thermoSys.nBoxes++;
			box[thermoSys.nBoxes-1].name = commands[1];
		}else if (commands[0] == "geometry") {box[thermoSys.nBoxes-1].geometry = commands[1];
		}else if (commands[0] == "surfacedensity") {box[thermoSys.nBoxes-1].solidDens = stod(commands[1]); // AA^-2
		}else if (commands[0] == "nlayersperwall") {box[thermoSys.nBoxes-1].nLayersPerWall = stoi(commands[1]);
		}else if (commands[0] == "distancebetweenlayers") {box[thermoSys.nBoxes-1].deltaLayers = stod(commands[1]);
		}else if (commands[0] == "size") {box[thermoSys.nBoxes-1].width[2] = stod(commands[1]); // AA. width[2] (z axis) always keeps for pore size.
		}else if (commands[0] == "fixvolume") {box[thermoSys.nBoxes-1].fix = true;} // AA. Used only for Gibbs ensemble.
		// Interaction parameters.
		else if (commands[2] == "vdwpotential"){
			ithBox = tls.FindIndex(box, thermoSys.nBoxes, commands[0]);
			ithSpecies = tls.FindIndex(fluid, thermoSys.nSpecies, commands[0]);
			jthSpecies = tls.FindIndex(fluid, thermoSys.nSpecies, commands[1]);
			if (ithSpecies >= 0 && jthSpecies >= 0) fluid[ithSpecies].vdwPot[jthSpecies] = commands[3];
			if (ithSpecies >= 0 && jthSpecies >= 0) fluid[jthSpecies].vdwPot[ithSpecies] = commands[3];
			if (ithBox >= 0 && jthSpecies >= 0) box[ithBox].fluid[jthSpecies].vdwPot[0] = commands[3];
		}else if (commands[2] == "sigma"){
			ithBox = tls.FindIndex(box, thermoSys.nBoxes, commands[0]);
			ithSpecies = tls.FindIndex(fluid, thermoSys.nSpecies, commands[0]);
			jthSpecies = tls.FindIndex(fluid, thermoSys.nSpecies, commands[1]);
			if (ithSpecies >= 0 && jthSpecies >= 0) fluid[ithSpecies].sigma[jthSpecies] = stod(commands[3]); // AA
			if (ithSpecies >= 0 && jthSpecies >= 0) fluid[jthSpecies].sigma[ithSpecies] = stod(commands[3]); // AA
			if (ithBox >= 0 && jthSpecies >= 0) box[ithBox].fluid[jthSpecies].sigma[0] = stod(commands[3]); // AA
		}else if (commands[2] == "epsilon"){
			ithBox = tls.FindIndex(box, thermoSys.nBoxes, commands[0]);
			ithSpecies = tls.FindIndex(fluid, thermoSys.nSpecies, commands[0]);
			jthSpecies = tls.FindIndex(fluid, thermoSys.nSpecies, commands[1]);
			if (ithSpecies >= 0 && jthSpecies >= 0) fluid[ithSpecies].epsilon[jthSpecies] = stod(commands[3]); // K
			if (ithSpecies >= 0 && jthSpecies >= 0) fluid[jthSpecies].epsilon[ithSpecies] = stod(commands[3]); // K
			if (ithBox >= 0 && jthSpecies >= 0) box[ithBox].fluid[jthSpecies].epsilon[0] = stod(commands[3]); // AA
		}else if (commands[2] == "rcut"){
			ithSpecies = tls.FindIndex(fluid, thermoSys.nSpecies, commands[0]);
			jthSpecies = tls.FindIndex(fluid, thermoSys.nSpecies, commands[1]);
			if (ithSpecies >= 0 && jthSpecies >= 0) fluid[ithSpecies].rcut[jthSpecies] = stod(commands[3]); // AA
			if (ithSpecies >= 0 && jthSpecies >= 0) fluid[jthSpecies].rcut[ithSpecies] = stod(commands[3]); // AA
		}else if (commands[2] == "numberofparticles"){
			ithBox = tls.FindIndex(box, thermoSys.nBoxes, commands[0]);
			jthSpecies = tls.FindIndex(fluid, thermoSys.nSpecies, commands[1]);
			if (ithBox >= 0 && jthSpecies >= 0) box[ithBox].fluid[jthSpecies].nParts = stoi(commands[3]);
		}
		// Sampling parameters.
		else if (commands[0] == "samplerdf"){
			ithSpecies = tls.FindIndex(fluid, thermoSys.nSpecies, commands[1]);
			jthSpecies = tls.FindIndex(fluid, thermoSys.nSpecies, commands[2]);
			sim.rdf[0] = ithSpecies;
			sim.rdf[1] = jthSpecies;
		}
		for (i=0; i<4; i++) commands[i] = "";
	}
	inFile.close();
	// Get remaining box widths and initial volumes according to the geometries given.
	maxRcut = 0.;
	for (i=0; i<thermoSys.nSpecies; i++){
		for (j=i; j<thermoSys.nSpecies; j++) if (maxRcut < fluid[i].rcut[j]) maxRcut = fluid[i].rcut[j];
	}
	for (i=0; i<thermoSys.nBoxes; i++){
		box[i].maxRcut = maxRcut;
		if (box[i].geometry == "sphere") {box[i].width[0] = box[i].width[1] = box[i].width[2];
		}else if (box[i].geometry == "cylinder"){
			box[i].width[0] = 2*maxRcut; // PBC along x axis.
			box[i].width[1] = box[i].width[2];
		}else if (box[i].geometry == "slit") {box[i].width[0] = box[i].width[1] = 2*maxRcut; // PBC along xy plane.
		}else box[i].width[0] = box[i].width[1] = box[i].width[2]; // Bulk phase.
		box[i].volume = ComputeVolume(box[i]);
		thermoSys.volume += box[i].volume;
	}
	// Get total number of particles and density per box.
	for (i=0; i<thermoSys.nBoxes; i++){
		for (j=0; j<thermoSys.nSpecies; j++) box[i].nParts += box[i].fluid[j].nParts;
		thermoSys.nParts += box[i].nParts;
	}
	// Set PBC restrictions.
	for (i=0; i<thermoSys.nBoxes; i++){
		if (box[i].geometry == "sphere"){
			box[i].PBC[0] = false;
			box[i].PBC[1] = false;
			box[i].PBC[2] = false;
		}else if (box[i].geometry == "cylinder"){
			box[i].PBC[0] = true;
			box[i].PBC[1] = false;
			box[i].PBC[2] = false;
		}else if (box[i].geometry == "slit"){
			box[i].PBC[0] = true;
			box[i].PBC[1] = true;
			box[i].PBC[2] = false;
		}else if (box[i].geometry == "bulk"){
			box[i].PBC[0] = true;
			box[i].PBC[1] = true;
			box[i].PBC[2] = true;
		} }
	// Set chemical potential to boxes if given.
	for (i=0; i<thermoSys.nBoxes; i++){
		for (j=0; j< thermoSys.nSpecies; j++) box[i].fluid[j].mu = fluid[j].mu;
	}
	// Set step size.
	//sim.dr = fluid[0].sigma[0]; //AA
}

