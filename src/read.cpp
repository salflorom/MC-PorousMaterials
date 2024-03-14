//Author: Santiago A. Flores Roman

#include <cmath> // cbrt, sqrt
#include <string> // string
#include <cstring> // strlen
#include <iostream> // cout, stod
#include <fstream> // ifstream
#include <cstdlib> // stod, stoi, strtol
#include <regex> // regex

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
long int MC::ReadInputFile(string inFileName){
	int i, j, ithBox, ithSpecies, jthSpecies;
	double maxRcut;
	string* commands;
	string line;
	Tools tls;
	ifstream inFile;
	long int set = 0;

	inFile.open(inFileName);
	if (!inFile.is_open()){ //Check if file was opened successfully.
		cout << "Error opening file: " << inFileName << endl;
		exit(EXIT_FAILURE);
	}
	while (getline(inFile, line)) {
		commands = tls.SplitString(line, ' ');
		for (i=0; i<4; i++) commands[i] = tls.LowerCase(commands[i]);
		// Simulation parameters.
		if (commands[0] == "projectname") sim.projName = commands[1];
		else if (commands[0] == "restart") sim.restart = true;
		else if (commands[0] == "productionsets") sim.nSets = stod(commands[1]);
		else if (commands[0] == "equilibriumsets") sim.nEquilSets = stod(commands[1]);
		else if (commands[0] == "stepsperset") sim.nStepsPerSet = stod(commands[1]);
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
		}else if (commands[0] == "printtrajectory") sim.printTrajectory = true;
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
	for (i=0; i<thermoSys.nBoxes; i++) sim.dr[i] = 0.1*box[i].width[2]; //AA
	// Set params according to read restart
	if (sim.restart) set = Restart();
	return set;
}
long int MC::Restart(void){
	char ch;
	int count;
	bool keepLooking;
	string lastLine, token;
	ifstream simFile;
	ostringstream inDirName, simFileName;
	size_t pos=0;
	int set=0;
	string substrings[11];
	string delimiter = "\t";

	inDirName << "./" << sim.projName;
	for (int i=0; i<thermoSys.nBoxes; i++){
		for (int j=0; j<thermoSys.nSpecies; j++){
			simFileName << inDirName.str() << "/" << box[i].name << "/simulation_" << fluid[j].name << ".log";
			simFile.open(simFileName.str());
			if (simFile.is_open()){ //Check if file was opened successfully.
				simFile.seekg(-2, ios_base::end); // go to one spot before the EOF
				keepLooking = true;
				while (keepLooking){
					simFile.get(ch); // Get current byte's data
					if ((int)simFile.tellg() <= 1){ // If the data was at or before the 0th byte
						simFile.seekg(0); // The first line is the last line
						keepLooking = false; // So stop there
					}else if (ch == '\n'){keepLooking = false; // Stop at the current position if the data was a newline.
					}else simFile.seekg(-2, ios_base::cur); // Move to the front of that data, then to the front of the data before it.
				}
				getline(simFile, lastLine); // Read the current line
				count = 0;
				while ((pos = lastLine.find(delimiter)) != string::npos){ // Split the line into an array of strings.
					token = lastLine.substr(0, pos);
					substrings[count] = token;
					lastLine.erase(0, pos + delimiter.length());
					count++;
				}
				if (substrings[0] == ""){
					cout << "Error reading log file. File: " << simFileName.str() << endl;
					exit(EXIT_FAILURE);
				}
				set = stol(substrings[0]);
				thermoSys.temp = stod(substrings[1]);
				box[i].width[2] = stod(substrings[2]);
				box[i].fluid[j].nParts = stoi(substrings[7]);
			}else{
				cout << "Warning: Last configuration not found. Path: " << simFileName.str() << endl;
				cout << "\tCurrent set, temperature, box width, and num. of particles will be set according to the input file." << endl << endl;
			}
			simFileName.str(string());
			simFileName.clear();
			simFile.close();
			simFile.clear();
		}
		box[i].nParts = 0;
		for (int j=0; j<thermoSys.nSpecies; j++) box[i].nParts += box[i].fluid[j].nParts;
	}
	thermoSys.nParts = 0;
	for (int i=0; i<thermoSys.nBoxes; i++) thermoSys.nParts += box[i].nParts;
	cout << "Restart simulation: Yes" << endl;
	cout << "Current system state: " << endl;
	cout << "\tCurrent set: " << set << endl;
	cout << "\tSystem temperature: " << thermoSys.temp << endl;
	for (int i=0; i<thermoSys.nBoxes; i++){
		cout << "\tBox " << box[i].name << " size: " << box[i].width[2] << endl;
		for (int j=0; j<thermoSys.nSpecies ; j++){
			cout << "\t\tNum. of particles for species " << fluid[j].name << " in the box: " << box[i].fluid[j].nParts << endl;
		}
		cout << "\tNum. of particles in the box: " << box[i].nParts << endl;
	}
	cout << endl;
	return set;
}

