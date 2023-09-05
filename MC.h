//Author: Santiago A. Flores Roman

#ifndef _MC_H_
#define _MC_H_

#include <string>
#include <cstdlib>
#include <cmath>
#include <random>

#define NBINS 200
#define MAXPART 9999 //Max. num. of particles. 0th part. is to save old config.
#define MAXSPECIES 2 //Max. num. of species.
#define MAXBOX 2 //Max. num. of boxes. 0th box is to save old comfiguration.
#define pi M_PI
#define kb 1.380649e-23 // J/K
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
			int acceptance, rejection, nDisplacements;
			int acceptanceVol, rejectionVol, nVolChanges;
			int acceptSwap, rejectSwap, nSwaps;
			int widomInsertions;
			double binWidth, widom[MAXBOX][MAXSPECIES], rdf[NBINS+1];
		} stats;
		struct Simulation{
			double nSets, nEquilSets, nCyclesPerSet, printEvery;
			double dr, dv;
			int nDispAttempts, nSwapAttempts, nVolAttempts;
			int cycle, rdf[2];
			string compute[3];
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
			bool PBC[3];
			double width[3];
			double solidDens, volume, deltaLayers;
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
		void AdjustMCMoves(void);
		void MoveParticle(void);
		void ChangeVolume(void);
		void RescaleCenterOfMass(Box, Box&);
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

		int* GetMCMoves(void);
		int GetNSets(void){return int(sim.nSets);}
		int GetNEquilSets(void){return int(sim.nEquilSets);}
		int GetNCyclesPerSet(void){return int(sim.nCyclesPerSet);}
		int GetPrintEvery(void){return int(sim.printEvery);}
		string* GetCompute(void){return sim.compute;}

		double NeighDistance(int, Particle, Particle);
		// Fluid-Fluid potentials //
		double HardSphere_Pot(int, int, int, int);
		double LJ_Pot(int, int, int, int);
		// EAM Ga potential vvvvv //
		double* EAMGA_Pot(int, int, int, int);
		double StepUnit(double, double, double);
		double EmbPot(double);
		double eDens(double);
		double PairPot(double);
		// EAM Ga potential ^^^^^ //
		// Fluid-Fluid potentials //

		// Solid-Fluid potentials //
		double SlitLJ(int, int, int);
		double CylindricalLJ(int, int, int);
		double SphericalLJ(int, int, int);
		// Solid-Fluid potentials //
/* **************************************************************************** */
};
/* **************************************************************************** */

#endif /* _MC_H_ */

