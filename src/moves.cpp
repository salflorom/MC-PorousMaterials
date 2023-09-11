//Author: Santiago A. Flores Roman

#include <cstring> // strlen
#include <iostream> // cout

#include "MC.h"
#include "tools.h"

using namespace std;

int* MC::GetMCMoves(void){
	static int moves[3];

	moves[0] = sim.nDispAttempts;
	moves[1] = sim.nVolAttempts;
	moves[2] = sim.nSwapAttempts;
	return moves;
}
void MC::InitialConfig(void){
	for (int i=0; i<thermoSys.nBoxes; i++){
		for (int j=0; j<thermoSys.nSpecies; j++){
			for (int k=1; k<=box[i].fluid[j].nParts; k++) InsertParticle(i, j, k);
		}
	}
	cout << "Initial configuration set." << endl;
	cout << endl;
}
void MC::InsertParticle(int ithBox, int ithSpecies, int index){
	double radius, theta, phi; // Spherical geometry.
	double rho; // Cylindrical geometry.

	if (box[ithBox].geometry == "sphere"){
		radius = Random()*box[ithBox].width[2]*0.5;
		theta = Random()*pi;
		phi = Random()*2*pi;
		box[ithBox].fluid[ithSpecies].particle[index].x = radius*sin(theta)*cos(phi);
		box[ithBox].fluid[ithSpecies].particle[index].y = radius*sin(theta)*sin(phi);
		box[ithBox].fluid[ithSpecies].particle[index].z = radius*cos(theta);
		// Frame is at (0,0,0). Moving part. to box center.
		box[ithBox].fluid[ithSpecies].particle[index].x += 0.5*box[ithBox].width[0];
		box[ithBox].fluid[ithSpecies].particle[index].y += 0.5*box[ithBox].width[1];
		box[ithBox].fluid[ithSpecies].particle[index].z += 0.5*box[ithBox].width[2];
	}else if (box[ithBox].geometry == "cylinder"){
		rho = Random()*box[ithBox].width[2]*0.5;
		phi = Random()*2*pi;
		box[ithBox].fluid[ithSpecies].particle[index].x = (Random()-0.5)*box[ithBox].width[0];
		box[ithBox].fluid[ithSpecies].particle[index].y = rho*sin(phi);
		box[ithBox].fluid[ithSpecies].particle[index].z = rho*cos(phi);
		// Frame is at (0,0,0). Moving part. to box center.
		box[ithBox].fluid[ithSpecies].particle[index].x += 0.5*box[ithBox].width[0];
		box[ithBox].fluid[ithSpecies].particle[index].y += 0.5*box[ithBox].width[1];
		box[ithBox].fluid[ithSpecies].particle[index].z += 0.5*box[ithBox].width[2];
	}else if (box[ithBox].geometry == "slit"){
		box[ithBox].fluid[ithSpecies].particle[index].x = Random()*box[ithBox].width[0];
		box[ithBox].fluid[ithSpecies].particle[index].y = Random()*box[ithBox].width[1];
		box[ithBox].fluid[ithSpecies].particle[index].z = Random()*box[ithBox].width[2];
	}else{
		box[ithBox].fluid[ithSpecies].particle[index].x = Random()*box[ithBox].width[0];
		box[ithBox].fluid[ithSpecies].particle[index].y = Random()*box[ithBox].width[1];
		box[ithBox].fluid[ithSpecies].particle[index].z = Random()*box[ithBox].width[2];
	}
}
void MC::AdjustMCMoves(void){
	double ratio;

	ratio = stats.acceptance*1./(1.*stats.nDisplacements);
	if (ratio < 0.2) sim.dr /= 1.1;
	else sim.dr *= 1.1;
	//if (sim.nVolAttempts > 0){
		//ratio = stats.acceptanceVol*1./(1.*stats.nVolChanges);
		//if (ratio < 0.3) sim.dv /= 1.1;
		//else sim.dv *= 1.1;
	//}
}
void MC::RescaleCenterOfMass(Box oldBox, Box& newBox){
	Tools tls;
	int i, j;
	double ratio = newBox.width[2]/oldBox.width[2];

	cout << ""; // To avoid core dump.
	for (i=0; i<thermoSys.nSpecies; i++){
		//if (pore.geometry == "sphere"){
			//for (i=1; i<=fluid.nParts; i++){
				//particle[i].x *= ratio;
				//particle[i].y *= ratio;
				//particle[i].z *= ratio;
			//}
		//}else if (pore.geometry == "cylinder"){
			//for (i=1; i<=fluid.nParts; i++){
				//particle[i].y *= ratio;
				//particle[i].z *= ratio;
			//}
		//}else if (pore.geometry == "slit"){
			//for (i=1; i<=fluid.nParts; i++) particle[i].z *= ratio;
		//}else{ //Bulk phase.
		for (j=1; j<=newBox.fluid[i].nParts; j++){
			newBox.fluid[i].particle[j].x *= ratio;
			newBox.fluid[i].particle[j].y *= ratio;
			newBox.fluid[i].particle[j].z *= ratio;
		}
	}
	//sim.dr *= ratio;
	//sim.dv *= ratio;
}
void MC::MoveParticle(void){
	Particle ithPart;
	int i, ithBox, ithSpecies, index, tmp;
	double rand;
	double oldEnergy, newEnergy, deltaEnergy, arg;
	double oldPairPotE, oldManyBodyE, newPairPotE, newManyBodyE;
	double oldBoxEnergy, newBoxEnergy;
	double deltaManyBodyE=0, deltaPairPotE=0, deltaBoxEnergy=0;

	stats.nDisplacements++;
	// Select particle.
	rand = int(Random()*thermoSys.nParts);
	tmp = 0;
	for (i=0; i<thermoSys.nBoxes; i++){
		tmp += box[i].nParts;
		if (rand <= tmp){
			ithBox = i;
			break;
		}
	}
	rand = int(Random()*box[ithBox].nParts);
	tmp = 0;
	for (i=0; i<thermoSys.nSpecies; i++){
		tmp += box[ithBox].fluid[i].nParts;
		if (rand <= tmp){
			ithSpecies = i;
			break;
		}
	}
	index = int(Random()*box[ithBox].fluid[ithSpecies].nParts)+1;
	ithPart = box[ithBox].fluid[ithSpecies].particle[index];
	box[ithBox].fluid[ithSpecies].particle[0] = ithPart; //Save old position of selected particle.
	EnergyOfParticle(ithBox, ithSpecies, index);
	oldManyBodyE = box[ithBox].fluid[ithSpecies].particle[index].manyBodyE;
	oldPairPotE = box[ithBox].fluid[ithSpecies].particle[index].pairPotE;
	oldBoxEnergy = box[ithBox].fluid[ithSpecies].particle[index].boxE;
	oldEnergy = box[ithBox].fluid[ithSpecies].particle[index].energy;
	// Move particle.
	ithPart.x += (Random()-0.5)*sim.dr; //AA
	ithPart.y += (Random()-0.5)*sim.dr; //AA
	ithPart.z += (Random()-0.5)*sim.dr; //AA
	PBC(ithBox, ithPart);
	box[ithBox].fluid[ithSpecies].particle[index] = ithPart;
	// Metropolis.
	EnergyOfParticle(ithBox, ithSpecies, index);
	newManyBodyE = box[ithBox].fluid[ithSpecies].particle[index].manyBodyE;
	newPairPotE = box[ithBox].fluid[ithSpecies].particle[index].pairPotE;
	newBoxEnergy = box[ithBox].fluid[ithSpecies].particle[index].boxE;
	newEnergy = box[ithBox].fluid[ithSpecies].particle[index].energy;
	deltaEnergy = newEnergy-oldEnergy;
	arg = deltaEnergy/thermoSys.temp;
	if (Random() < exp(-arg)){
		stats.acceptance++;
		deltaManyBodyE = newManyBodyE - oldManyBodyE;
		deltaPairPotE = newPairPotE - oldPairPotE;
		deltaBoxEnergy = newBoxEnergy - oldBoxEnergy;
		box[ithBox].manyBodyE += deltaManyBodyE;
		box[ithBox].pairPotE += 0.5*deltaPairPotE;
		box[ithBox].boxE += deltaBoxEnergy;
		box[ithBox].energy += deltaManyBodyE + 0.5*deltaPairPotE + deltaBoxEnergy;
	}else{
		stats.rejection++;
		ithPart = box[ithBox].fluid[ithSpecies].particle[0];
		box[ithBox].fluid[ithSpecies].particle[index] = ithPart;
	}
}
// Performs MC trial move: change of volume.
// Note: It applies either for Gibbs or NPT ensemble.
// Note 2: For NPT ensemble, it only works in box 1.
// Note 3: For Gibbs ensemble, it only works with 2 boxes.
void MC::ChangeVolume(void){
	Box oldBox0, oldBox1;
	int ithBox;
	double deltaE0, deltaE1;
	double logNewVol, deltaVol0, deltaVol1;
	double arg0, arg1;
	double extPress = thermoSys.press/kb*1e-30; // K/AA^3

	stats.nVolChanges++;
	if (extPress > 0){ // NPT or Gibbs-NPT ensemble.
		// Record current config.
		do{
			ithBox = round(Random()*thermoSys.nBoxes);
		}while (box[ithBox].fix);
		oldBox0 = box[ithBox]; // Record old information.
		logNewVol = log(box[ithBox].volume) + (Random()-0.5)*sim.dv; //Perform volume change step.
		box[ithBox].volume = exp(logNewVol);
		RescaleCenterOfMass(oldBox0, box[ithBox]); //Rescale to trial config.
		// Set box sizes according to trial volume.
		//if (pore.geometry == "sphere") pore.width[0] = pore.width[1] = pore.width[2] = finalBoxWidth;
		//else if (pore.geometry == "cylinder") pore.width[1] = pore.width[2] = finalBoxWidth;
		//else if (pore.geometry == "slit") pore.width[2] = finalBoxWidth;
		//else pore.width[0] = pore.width[1] = pore.width[2] = finalBoxWidth;
		box[ithBox].width[2] = ComputeBoxWidth(box[ithBox], box[ithBox].volume);
		box[ithBox].width[0] = box[ithBox].width[1] = box[ithBox].width[2];
		// Compute trial energy.
		BoxEnergy(ithBox);
		deltaE0 = box[ithBox].energy-oldBox0.energy;
		deltaVol0 = box[ithBox].volume-oldBox0.volume;
		arg0 = deltaE0 + extPress*deltaVol0 - (box[ithBox].nParts+1)*log(box[ithBox].volume/oldBox0.volume)*thermoSys.temp;
		if (Random() > exp(-arg0/thermoSys.temp)){ //Rejected trial move.
			stats.rejectionVol++;
			box[ithBox] = oldBox0;
		}else{ // Accepted trial move. Keep trial config.
			stats.acceptanceVol++;
			thermoSys.volume += deltaVol0;
		}
	}else{ // Gibbs-NVT ensemble.
		oldBox0 = box[0]; // Record old information.
		oldBox1 = box[1];
		logNewVol = log(box[0].volume/box[1].volume) + (Random()-0.5)*sim.dv; //Perform volume change step.
		box[0].volume = thermoSys.volume*exp(logNewVol)/(1+exp(logNewVol));
		box[1].volume = thermoSys.volume-box[0].volume;
		RescaleCenterOfMass(oldBox0, box[0]); //Rescale to trial config.
		RescaleCenterOfMass(oldBox1, box[1]); //Rescale to trial config.
		// Set box sizes according to trial volume.
		box[0].width[2] = ComputeBoxWidth(box[0], box[0].volume);
		box[0].width[0] = box[0].width[1] = box[0].width[2];
		box[1].width[2] = ComputeBoxWidth(box[1], box[1].volume);
		box[1].width[0] = box[1].width[1] = box[1].width[2];
		// Compute trial energy.
		BoxEnergy(0);
		BoxEnergy(1);
		deltaE0 = box[0].energy-oldBox0.energy;
		deltaE1 = box[1].energy-oldBox1.energy;
		arg0 = deltaE0 + (box[0].nParts+1)*log(box[0].volume/oldBox0.volume)*thermoSys.temp;
		arg1 = deltaE1 + (box[1].nParts+1)*log(box[1].volume/oldBox1.volume)*thermoSys.temp;
		if (Random() > exp(-(arg0+arg1)/thermoSys.temp)){ //Rejected trial move.
			stats.rejectionVol++;
			box[0] = oldBox0;
			box[1] = oldBox1;
		}else stats.acceptanceVol++; // Accepted trial move. Keep trial config.
	}
	for (int i=0; i<thermoSys.nBoxes; i++){
		if (box[i].geometry == "bulk" && box[i].width[2] < box[i].maxRcut){
			cout << "Warning: Box surpassed the minimum allowed size: ";
			cout << box[i].maxRcut << endl;
		}
	}
}
// Performs MC trial move: change of volume.
// Note: It applies either for Gibbs or muVT ensemble.
// Note 2: For muVT ensemble, it only works in box 1.
// Note 3: For Gibbs ensemble, it only works with 2 boxes.
void MC::SwapParticle(void){
	Tools tls;
	int tmp, ithSpecies, inIndex, outIndex, inBox, outBox;
	double particleMass, thermalWL, zz;
	double inEnergy, outEnergy, deltaE, pairPotE, manyBodyE, boxE;
	double rand;
	double arg=0, temp=thermoSys.temp;

	if (thermoSys.nBoxes == 1){
		rand = int(Random()*box[0].nParts);
		tmp = 0;
		for (int i=0; i<thermoSys.nSpecies; i++){
			tmp += box[0].fluid[i].nParts;
			if (rand <= tmp){
				ithSpecies = i;
				break;
			}
		}
		particleMass = fluid[ithSpecies].molarMass/na*1e-3; // kg/particle
		thermalWL = planck/sqrt(2.0*pi*particleMass*kb*thermoSys.temp)*1e10; // AA
		zz = exp(fluid[ithSpecies].mu/thermoSys.temp)/tls.Pow(thermalWL,3); // AA^-3
		if (Random() < 0.5) { // Try inserting particle.
			stats.nSwaps++;
			// Insert particle at random position.
			inIndex = box[0].fluid[ithSpecies].nParts+1;
			InsertParticle(0, ithSpecies, inIndex);
			EnergyOfParticle(0, ithSpecies, inIndex);
			inEnergy = box[0].fluid[ithSpecies].particle[inIndex].energy; // K
			arg = zz*box[0].volume*exp(-inEnergy/thermoSys.temp)/(box[0].nParts+1); // Acceptance criterion (for inserting particle).
			if (Random() < arg){
				stats.acceptSwap++; // Accepted: Insert particle.
				thermoSys.nParts++;
				box[0].nParts++;
				box[0].fluid[ithSpecies].nParts++;
				manyBodyE = box[0].fluid[ithSpecies].particle[inIndex].manyBodyE;
				pairPotE = 0.5*box[0].fluid[ithSpecies].particle[inIndex].pairPotE;
				boxE = box[0].fluid[ithSpecies].particle[inIndex].boxE;
				box[0].manyBodyE += manyBodyE; // K
				box[0].pairPotE += pairPotE; // K
				box[0].boxE += boxE; // K
				box[0].energy += manyBodyE + pairPotE + boxE;
			}else stats.rejectSwap++;
		}else if (box[0].nParts > 0){ // Try removing particle (only if there are particles in the box).
			stats.nSwaps++;
			// Select random particle.
			outIndex = int(Random()*(box[0].fluid[ithSpecies].nParts)+1);
			EnergyOfParticle(0, ithSpecies, outIndex); //K
			outEnergy = box[0].fluid[ithSpecies].particle[outIndex].energy; // K
			arg = box[0].nParts*exp(outEnergy/thermoSys.temp)/(zz*box[0].volume); // Acceptance criterion (for removing particle).
			if (Random() < arg){ // Accepted: Remove particle.
				stats.acceptSwap++;
				box[0].fluid[ithSpecies].particle[outIndex] = box[0].fluid[ithSpecies].particle[box[0].fluid[ithSpecies].nParts];
				thermoSys.nParts--;
				box[0].nParts--;
				box[0].fluid[ithSpecies].nParts--;
				manyBodyE = box[0].fluid[ithSpecies].particle[outIndex].manyBodyE;
				pairPotE = 0.5*box[0].fluid[ithSpecies].particle[outIndex].pairPotE;
				boxE = box[0].fluid[ithSpecies].particle[outIndex].boxE;
				box[0].manyBodyE -= manyBodyE; // K
				box[0].pairPotE -= pairPotE; // K
				box[0].boxE -= boxE; // K
				box[0].energy -= manyBodyE + pairPotE + boxE;
			}else stats.rejectSwap++;
		}
	}else{ // Gibbs ensemble.
		stats.nSwaps++;
		if (Random() < 0.5){
			inBox = 0;
			outBox = 1;
		}else{
			inBox = 1;
			outBox = 0;
		}
		if (box[outBox].nParts > 0){
			rand = int(Random()*box[inBox].nParts);
			tmp = 0;
			for (int i=0; i<thermoSys.nSpecies; i++){
				tmp += box[inBox].fluid[i].nParts;
				if (rand <= tmp){
					ithSpecies = i;
					break;
				}
			}
			inIndex = box[inBox].fluid[ithSpecies].nParts+1;
			InsertParticle(inBox, ithSpecies, inIndex);
			EnergyOfParticle(inBox, ithSpecies, inIndex);
			inEnergy = box[inBox].fluid[ithSpecies].particle[inIndex].energy; // K
			stats.widom[inBox][ithSpecies] += box[inBox].volume*exp(-inEnergy/thermoSys.temp)/(box[inBox].nParts+1); // Compute mu through Widom.
			stats.widomInsertions++;
			outIndex = int(Random()*box[outBox].fluid[ithSpecies].nParts + 1);
			EnergyOfParticle(outBox, ithSpecies, outIndex); //K
			outEnergy = box[outBox].fluid[ithSpecies].particle[outIndex].energy; // K
			deltaE = inEnergy-outEnergy;
			arg = deltaE+log(box[outBox].volume*(box[inBox].nParts+1)/(box[inBox].volume*box[outBox].nParts))*temp; // Acceptance criterion (for removing particle).
			if (Random() < exp(-arg/temp)){
				stats.acceptSwap++; // Accepted: Insert particle.
				box[inBox].nParts++;
				box[inBox].fluid[ithSpecies].nParts++;
				manyBodyE = box[inBox].fluid[ithSpecies].particle[inIndex].manyBodyE;
				pairPotE = 0.5*box[inBox].fluid[ithSpecies].particle[inIndex].pairPotE;
				boxE = box[inBox].fluid[ithSpecies].particle[inIndex].boxE;
				box[inBox].manyBodyE += manyBodyE; // K
				box[inBox].pairPotE += pairPotE; // K
				box[inBox].boxE += boxE; // K
				box[inBox].energy += manyBodyE + pairPotE + boxE;
				box[outBox].fluid[ithSpecies].particle[outIndex] = box[outBox].fluid[ithSpecies].particle[box[outBox].fluid[ithSpecies].nParts];
				box[outBox].nParts--;
				box[outBox].fluid[ithSpecies].nParts--;
				manyBodyE = box[outBox].fluid[ithSpecies].particle[outIndex].manyBodyE;
				pairPotE = 0.5*box[outBox].fluid[ithSpecies].particle[outIndex].pairPotE;
				boxE = box[outBox].fluid[ithSpecies].particle[outIndex].boxE;
				box[outBox].manyBodyE -= manyBodyE; // K
				box[outBox].pairPotE -= pairPotE; // K
				box[outBox].boxE -= boxE; // K
				box[outBox].energy -= manyBodyE + pairPotE + boxE;
			}else stats.rejectSwap++;
		}
	}
}
// Performs Widom insertions to compute mu in each box.
// Note: It only works for simple fluids.
// Note 2: The Widom parameter is computed randomly in any box unless the chosen box is emtpy.
void MC::ComputeWidom(void){
	Tools tls;
	double volume;
	double energy=0.;

	if (thermoSys.nBoxes == 1 && sim.nSwapAttempts == 0){
		InsertParticle(0, 0, MAXPART-1); // Insert virtual particle.
		stats.widomInsertions++;
		EnergyOfParticle(0, 0, MAXPART-1);
		energy = box[0].fluid[0].particle[MAXPART-1].energy;
		if (sim.nVolAttempts > 0){
			stats.widom[0][0] += box[0].volume*exp(-energy/thermoSys.temp)/(box[0].fluid[0].nParts+1);
		}else stats.widom[0][0] += exp(-energy/thermoSys.temp);
		stats.widom[0][0] += exp(-energy/thermoSys.temp);
	}
}

