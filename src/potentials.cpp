//Author: Santiago A. Flores Roman

#include <cstring> // strlen
#include <cmath> // pi, exp
#include <climits> // INT_MAX

#include "MC.h"
#include "tools.h"
#include "NumRecipes/hypgeo.h" // hypgeo

using namespace std;


double MC::NeighDistance(int ithBox, Particle ithPart, Particle jthPart){
	double dist, dx, dy, dz;

	dx = ithPart.x - jthPart.x; //AA
	dy = ithPart.y - jthPart.y; //AA
	dz = ithPart.z - jthPart.z; //AA
	if (box[ithBox].PBC[0]) dx -= round(dx/box[ithBox].width[0]) * box[ithBox].width[0]; //AA
	if (box[ithBox].PBC[1]) dy -= round(dy/box[ithBox].width[1]) * box[ithBox].width[1]; //AA
	if (box[ithBox].PBC[2]) dz -= round(dz/box[ithBox].width[2]) * box[ithBox].width[2]; //AA
	dist = sqrt(dx*dx + dy*dy + dz*dz); //AA
	return dist; //AA
}

// ---------- Fluid-Fluid potentials ---------- //
double MC::HardSphere_Pot(int ithBox, int ithSpecies, int jthSpecies, int index){
	Tools tls;
	Particle ithPart, jthPart;
	double rij, sig;
	double energy=0;

	ithPart = box[ithBox].fluid[ithSpecies].particle[index];
	sig = fluid[ithSpecies].sigma[jthSpecies]; // AA
	for (int j=1; j<=box[ithBox].fluid[jthSpecies].nParts; j++){
		jthPart = box[ithBox].fluid[jthSpecies].particle[j];
		rij = NeighDistance(ithBox, ithPart, jthPart);
		if (rij <= sig) energy += INT_MAX;
	}
	return energy; //K
}
double MC::LJ_Pot(int ithBox, int ithSpecies, int jthSpecies, int index){
	Tools tls;
	Particle ithPart, jthPart;
	double rij, eps, sig, rcut;
	double energy=0;

	ithPart = box[ithBox].fluid[ithSpecies].particle[index];
	eps = fluid[ithSpecies].epsilon[jthSpecies]; // K
	sig = fluid[ithSpecies].sigma[jthSpecies]; // AA
	rcut = fluid[ithSpecies].rcut[jthSpecies]; // AA
	for (int j=1; j<=box[ithBox].fluid[jthSpecies].nParts; j++){
		jthPart = box[ithBox].fluid[jthSpecies].particle[j];
		if (ithPart != jthPart){
			rij = NeighDistance(ithBox, ithPart, jthPart);
			if (rij <= rcut) energy += 4.*eps*(tls.Pow(sig/rij,12)-tls.Pow(sig/rij,6));
		}
	}
	return energy; //K
}
// EAM Ga potential vvv //
// Paper:
// Belashchenko, D.K., 2012.
// Computer Simulation of the Properties of Liquid Metals: Gallium, Lead, and Bismuth.
// Russ. J. Phys. Chem. A, 86, pp.779-790.
double* MC::EAMGA_Pot(int ithBox, int ithSpecies, int jthSpecies, int index){
	Particle ithPart, jthPart;
	static double energy[2];
	double rij, rcut;
	double rho=0, phiLow=0;
	double eVToK = (801088317./5.0e27)/kb;

	ithPart = box[ithBox].fluid[ithSpecies].particle[index];
	rcut = fluid[ithSpecies].rcut[jthSpecies]; // AA
	for (int j=1; j<=box[ithBox].fluid[jthSpecies].nParts; j++){
		jthPart = box[ithBox].fluid[jthSpecies].particle[j];
		if (ithPart != jthPart){
			rij = NeighDistance(ithBox, ithPart, jthPart);
			if (rij <= rcut){
				rho += eDens(rij);
				phiLow += PairPot(rij);
			}
		}
	}
	energy[0] = EmbPot(rho)*eVToK;
	energy[1] = phiLow*eVToK;
	return energy;
}
double MC::StepUnit(double radius, double leftLim, double rightLim){
	if (leftLim <= radius && radius < rightLim) return 1.;
	else return 0.;
}
double MC::EmbPot(double rho){
	double rhoIntervals[7] = {1.00000,  0.92000,  0.870000,  0.800000,  0.750000,  0.650000,  1.400000};
	double aValues[7]      = {0.00000, -1.91235, -1.904030, -1.897380, -1.883520, -1.852620, -1.822820};
	double bValues[7]      = {0.00000,  0.00000, -0.208000, -0.058000, -0.338000, -0.898000,  0.302000};
	double cValues[7]      = {0.00000,  1.30000, -1.500000,  2.000000,  5.600000, -6.000000,  2.000000};
	double phi = 0;

	if (rhoIntervals[1] <= rho && rho <= rhoIntervals[6]) phi = aValues[1] + cValues[1]*(rho-rhoIntervals[0])*(rho-rhoIntervals[0]);
	for (int i=2; i<6; i++){
		if (rhoIntervals[i] <= rho && rho <= rhoIntervals[i-1]){
			phi = aValues[i] + bValues[i]*(rho-rhoIntervals[i-1]) + cValues[i]*(rho-rhoIntervals[i-1])*(rho-rhoIntervals[i-1]);
			break;
		}
	}
	if (rho <= rhoIntervals[5]){
		phi = (aValues[6] + bValues[6]*(rho-rhoIntervals[5]) + cValues[6]*(rho-rhoIntervals[5])*(rho-rhoIntervals[5])) * (2*rho/rhoIntervals[5]-(rho/rhoIntervals[5])*(rho/rhoIntervals[5]));
	}
	return phi;
}
double MC::eDens(double radius){
	double pValues[3] = {0, 2.24450, 1.2};

	return pValues[1]*exp(-pValues[2]*radius);
}
double MC::PairPot(double radius){
	double rIntervals[7] = {0.0, 2.15, 2.75, 3.35, 4.00, 6.50, 8.30};
	double aValues[9][6] = {{0.0, -0.65052509307861e-01, -0.15576396882534e+00, -0.13794735074043e+00, 0.13303710147738e-01,  0.00000000000000e+00},
							{0.0, -0.32728102803230e+00, -0.16365580260754e-01,  0.78778542578220e-01, 0.59769893996418e-02,  0.00000000000000e+00},
							{0.0,  0.51590444127493e+01,  0.20955204046244e+00, -0.83622260891495e-01, 0.57411338894840e-01, -0.60454444423660e-02},
							{0.0,  0.90195221829217e+02, -0.97550604734748e+00, -0.44410858010987e+01, 0.19517888219051e+00, -0.13258585494287e+00},
							{0.0,  0.72322004859499e+03, -0.11625479189815e+02, -0.36415106938231e+02, 0.32162310059276e+00, -0.34988482891053e+00},
							{0.0,  0.27788989409594e+04, -0.58549935696765e+02, -0.13414583419234e+03, 0.30195698240893e+00, -0.45183606796559e+00},
							{0.0,  0.56037895713613e+04, -0.15186293377510e+03, -0.25239146992011e+03, 0.14850603977640e+00, -0.31733856650298e+00},
							{0.0,  0.57428084950480e+04, -0.19622924502226e+03, -0.23858760191913e+03, 0.36233874262589e-01, -0.11493645479281e+00},
							{0.0,  0.23685488320885e+04, -0.98789413798382e+02, -0.90270667293646e+02, 0.34984220138018e-02, -0.16768950999376e-01}};
	double phi=0;

	if (rIntervals[1] < radius && radius <= rIntervals[6]){
		for (int i=1; i<6; i++){
			for (int m=0; m<9; m++){
				phi += aValues[m][i] * pow((radius-rIntervals[i+1]),m) * StepUnit(radius, rIntervals[i], rIntervals[i+1]);
			}
		}
	}else if (rIntervals[0] < radius && radius <= rIntervals[1]){
		phi = 0.619588 - 51.86268*(2.15-radius) + 27.8*(exp(1.96*(2.15-radius))-1);
	}
	return phi;
}
// EAM Ga potential ^^^ //
// ---------- Fluid-Fluid potentials ---------- //

// ---------- Solid-Fluid potentials ---------- //
// Paper:
// Steele, W.A., 1973.
// The Physical Interaction of Gases with Crystalline Solids: I. Gas-Solid Energies and Properties of Isolated Adsorbed Atoms.
// Surf. Sci., 36(1), pp.317-352.
// Steele 10-4-3 potential extended by Jason for multiple layers.
// Domain of potential: z in [0,H], where H is the pore size.
double MC::SlitLJ_Pot(int ithBox, int ithSpecies, int index){
	Tools tls;
	double z, t1, t2, t3, t4;
	int nLayers = box[ithBox].nLayersPerWall;
	double sizeR = box[ithBox].width[2]*0.5; // AA
	double dens = box[ithBox].solidDens; // AA^-2
	double eps = box[ithBox].fluid[ithSpecies].epsilon[0]; // K
	double sig = box[ithBox].fluid[ithSpecies].sigma[0]; // AA
	double delta = box[ithBox].deltaLayers; // AA
	double usf = 0;

	z = box[ithBox].fluid[ithSpecies].particle[index].z-sizeR;
	for (int i=0; i<nLayers; i++){
		t1 = tls.Pow(sig/(sizeR+i*delta+z),10);
		t2 = tls.Pow(sig/(sizeR+i*delta-z),10);
		t3 = tls.Pow(sig/(sizeR+i*delta+z),4);
		t4 = tls.Pow(sig/(sizeR+i*delta-z),4);
		usf += 0.2*(t1+t2)-0.5*(t3+t4);
	}
	usf *= 4.*pi*eps*dens*sig*sig; // K
	return usf;
}
// Paper:
// Tjatjopoulos et al., 1988.
// Molecule-Micropore Interaction Potentials.
// J. Phys. Chem., 92(13), pp.4006-4007.
// Source for hypergeometric function:
// Press, W.H., 2007.
// Numerical recipes 3rd edition: The art of scientific computing.
// Cambridge university press.
double MC::CylindricalLJ_Pot(int ithBox, int ithSpecies, int index){
	Tools tls;
	double yPos, zPos;
	double x, u1, u2, usf, rho;
	double sizeR = box[ithBox].width[2]*0.5; // AA
	double dens = box[ithBox].solidDens; // AA^-2
	double eps = box[ithBox].fluid[ithSpecies].epsilon[0]; // K
	double sig = box[ithBox].fluid[ithSpecies].sigma[0]; // AA

	yPos = box[ithBox].fluid[ithSpecies].particle[index].y - 0.5*box[ithBox].width[1];
	zPos = box[ithBox].fluid[ithSpecies].particle[index].z - 0.5*box[ithBox].width[2];
	rho = sqrt(tls.Pow(yPos,2) + tls.Pow(zPos,2));
	// Reducing units/
	dens *= sig*sig;
	sizeR /= sig;
	rho /= sig;
	x = sizeR-rho;
	u1 = 63. /32. * 1./(tls.Pow(x,10)*tls.Pow(2.0-x/sizeR,10)) * hypgeo(-4.5,-4.5,1.,tls.Pow(1-x/sizeR,2));
	u2 = 3. / (tls.Pow(x,4)*tls.Pow(2.0-x/sizeR,4)) * hypgeo(-1.5,-1.5,1.,tls.Pow(1.0-x/sizeR,2));
	usf = (u1-u2)*pi*pi*dens*eps; // K
	return usf;
}
// Paper:
// Baksh, M.S.A. and Yang, R.T., 1991.
// Model for Spherical Cavity Radii and Potential Functions of Sorbates in Zeolites.
// AIChE J., 37(6), pp.923-930.
double MC::SphericalLJ_Pot(int ithBox, int ithSpecies, int index){
	Tools tls;
	double xPos, yPos, zPos;
	double r, x, usf;
	double sizeR = box[ithBox].width[2]*0.5;
	double dens = box[ithBox].solidDens;
	double eps = box[ithBox].fluid[ithSpecies].epsilon[0];
	double sig = box[ithBox].fluid[ithSpecies].sigma[0];
	double u1=0, u2=0;

	xPos = box[ithBox].fluid[ithSpecies].particle[index].x - 0.5*box[ithBox].width[0];
	yPos = box[ithBox].fluid[ithSpecies].particle[index].y - 0.5*box[ithBox].width[1];
	zPos = box[ithBox].fluid[ithSpecies].particle[index].z - 0.5*box[ithBox].width[2];
	r = sqrt(tls.Pow(xPos,2) + tls.Pow(yPos,2) + tls.Pow(zPos,2));
	// Reducing units/
	dens *= sig*sig;
	sizeR /= sig;
	r /= sig;
	x = sizeR-r;
	for (int i=0; i<10; i++){
		u1 += 1. / tls.Pow(sizeR,i) / tls.Pow(x,10-i);
		u1 += tls.Pow(-1.0,i) / tls.Pow(sizeR,i) / tls.Pow(x-2.0*sizeR,10-i);
		if (i < 4){
			u2 += 1. / tls.Pow(sizeR,i) / tls.Pow(x,4-i);
			u2 += tls.Pow(-1.0,i) / tls.Pow(sizeR,i) / tls.Pow(x-2.0*sizeR,4-i);
		}
	}
	usf = (2./5.) * u1 - u2;
	usf *= 2.*pi*eps*dens; //K
	return usf;
}
// ---------- Solid-Fluid potentials ---------- //

