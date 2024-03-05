//Author: Santiago A. Flores Roman

#ifndef _MC_H_
#define _MC_H_

#include <string> // string
#include <cmath> //M_PI
#include <random> // random_device, mt19937_64, uniform_real_distribution

#define NBINS 100
#define MAXPART 9999 //Max. num. of particles. 0th part. is to save old config.
#define MAXSPECIES 2 //Max. num. of species.
#define MAXBOX 2 //Max. num. of boxes. 0th box is to save old comfiguration.
#define pi M_PI
#define kb  1.380649e-23 // J/K
#define na 6.02214076e23 // mol^-1
#define planck 6.62607015e-34 // J*s

using namespace std;

/* **************************************************************************** */
class MC {
/* **************************************************************************** */
	protected:
/* **************************************************************************** */
		// Use random_device to generate a seed for Mersenne twister engine.
		random_device rd{};
		// Use Mersenne twister engine to generate pseudo-random numbers.
		mt19937_64 engine{rd()};
		// "Filter" MT engine's output to generate pseudo-random double values,
		// **uniformly distributed** on the interval [0, 1).
		uniform_real_distribution<double> dis{0.0,0.9999};

		struct Stats{
			int acceptance[MAXBOX], rejection[MAXBOX], nDisplacements[MAXBOX];
			int acceptanceVol[MAXBOX], rejectionVol[MAXBOX], nVolChanges;
			int acceptSwap, rejectSwap, nSwaps;
			int widomInsertions[MAXBOX];
			double binWidth, widom[MAXBOX][MAXSPECIES], rdf[NBINS+1];
		} stats;
		struct Simulation{
			bool printTrajectory;
			double nSets, nEquilSets, nStepsPerSet, printEvery;
			double dr[MAXBOX], dv[MAXBOX];
			int nDispAttempts, nSwapAttempts, nVolAttempts;
			int cycle, rdf[2];
			string projName;
		} sim;
		struct Particle{
			double x, y, z;
			double energy, manyBodyE, pairPotE, boxE; // //Particle energy.
			bool operator!=(Particle bPart){
				if ((this->x != bPart.x) || (this->y != bPart.y) || (this->z != bPart.z)){
					return true;
				}else return false;
			}
		};
		struct Fluid{
			string name;
			string vdwPot[MAXSPECIES];
			int nParts;
			double molarMass, mu, muEx;
			double sigma[MAXSPECIES], epsilon[MAXSPECIES], rcut[MAXSPECIES];
			//Energy of the ith species in the jth box (ignoring interactions with the box).
			Particle particle[MAXPART];
		} fluid[MAXSPECIES]; //For single species properties.
		struct Box{
			string name, geometry, vdwPot; //pot: solid-fluid potential.
			int nParts;
			bool fix, PBC[3];
			double width[3];
			double oldEnergy, solidDens, volume, deltaLayers, maxRcut;
			double energy, manyBodyE, pairPotE, boxE; // Energy of the box.
			int nLayersPerWall;
			Fluid fluid[MAXSPECIES]; //For solid-species properties.
		} box[MAXBOX];
		struct ThermodynamicSystem{
			double temp, volume, press;
			int nBoxes, nSpecies, nParts;
		} thermoSys;
/* **************************************************************************** */
	public:
/* **************************************************************************** */
		MC(void);
		void OutputDirectory(void);
		void ReadInputFile(string);
		void ResetStats(void);
		void PrintParams(void);
		void InsertParticle(int, int, int);
		void InitialConfig(void);
		void MinimizeEnergy(void);
		void AdjustMCMoves(int);
		void MoveParticle(void);
		void ChangeVolume(void);
		void RescaleCenterOfMass(Box, Box&, int);
		void SwapParticle(void);
		void PBC(int, Particle&);
		void ComputeRDF(void);
		void PrintRDF(void);
		void Metropolis(int);
		void PrintTrajectory(int);
		void PrintLog(int);
		void PrintStats(int);
		void ComputeWidom(void);
		void ComputeChemicalPotential(void);
		void EnergyOfParticle(int, int, int);
		double Random(void){return dis(engine);} //random num. in the interval [0,1).
		double ComputeVolume(Box);
		double ComputeBoxWidth(Box, double);
		void BoxEnergy(int);
		void CorrectEnergy(void);

		int* GetMCMoves(void);
		int GetNSets(void){return int(sim.nSets);}
		int GetNEquilSets(void){return int(sim.nEquilSets);}
		int GetNStepsPerSet(void){return int(sim.nStepsPerSet);}
		int GetPrintEvery(void){return int(sim.printEvery);}

		double NeighDistance(int, Particle, Particle);
		double StepUnit(double, double, double);
		// Fluid-Fluid potentials //
		double HardSphere_Pot(int, int, int, int);
		// Lennard-Jones 12-6 potential //
		double LJ126_Pot(int, int, int, int);
		// EAM Na potential //
		double* EAMNa_Pot(int, int, int, int);
		double EAMNa_EmbPot(double);
		double EAMNa_eDens(double);
		double EAMNa_PairPot(double);
		// EAM K potential //
		double* EAMK_Pot(int, int, int, int);
		double EAMK_EmbPot(double);
		double EAMK_eDens(double);
		double EAMK_PairPot(double);
		// EAM Rb potential //
		double* EAMRb_Pot(int, int, int, int);
		double EAMRb_EmbPot(double);
		double EAMRb_eDens(double);
		double EAMRb_PairPot(double);
		// EAM Ga potential //
		double* EAMGa_Pot(int, int, int, int);
		double EAMGa_EmbPot(double);
		double EAMGa_eDens(double);
		double EAMGa_PairPot(double);
		// Fluid-Fluid potentials //

		// Solid-Fluid potentials //
		double SlitLJ_Pot(int, int, int);
		double CylindricalLJ10_4(int, int, int);
		double CylindricalSteele10_4_3(int, int, int);
		double SphericalLJ_Pot(int, int, int);
		// Solid-Fluid potentials //
/* **************************************************************************** */
};
/* **************************************************************************** */

#endif /* _MC_H_ */

