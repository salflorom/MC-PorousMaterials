//Author: Santiago A. Flores Roman

#ifndef _MC_H_
#define _MC_H_

#include <string>
#include <cstdlib>
#include <cmath>
#include <vector>
#include <random>

#define NBINS 100
#define MAXPART 99999 //Max. num. of particles.
#define pi M_PI
#define kb 1.380649e-23 // J/K
#define na 6.02214e23 // mol^-1
#define plank 1.054571817e-34 // J*s

using namespace std;

/* **************************************************************************** */
class MC {
/* **************************************************************************** */
	public:
/* **************************************************************************** */
		MC(void);
		void ReadInputFile(string inFileName);
		void ResetStats(void);
		void PrintParams(void);
		void ComputeBoxSize(void);
		void InitialConfig(void);
		void MinimizeEnergy(void);
		void MoveParticle(void);
		void ChangeVolume(void);
		void ExchangeParticle(void);
		void PBC(int index);
		void ComputeRDF(void);
		void PrintRDF(int set);
		void Metropolis(int index);
		void ResetParticle(int index);
		void CreateEXYZ(int step);
		void CreateLogFile(int step);
		void PrintStats(int step);
		void ComputeWidom(void);
		void ComputeChemicalPotential(void);
		double SystemEnergy(void);
		double EnergyOfParticle(int index);
		double LJ_Energy(int i, int j);
		// For EAM Ga potential vvvvv //
		double EAMGA_Energy(int index);
		double StepUnit(double radius, double leftLim, double rightLim);
		double EmbPot(double rho);
		double eDens(double radius);
		double PairPot(double radius);
		// For EAM Ga potential ^^^^^ //
		double NeighDistance(int i, int j);
		double Random(void){return dis(engine);} //random num. in the interval [0,1).
		double* GetMCMoveProbabilities(void);
		int GetNSets(void){return int(sim.nSets);}
		int GetNEquilSets(void){return int(sim.nEquilSets);}
		int GetNStepsPerSet(void){return int(sim.nStepsPerSet);}
		int GetPrintEvery(void){return int(sim.printEvery);}
		string* GetCompute(void){return sim.compute;}
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
			int acceptInsertion, rejectInsertion;
			int acceptDeletion, rejectDeletion, nExchanges;
			int widomInsertions;
			double binWidth, widom, rdf[NBINS+1];
		} stats;
		struct Simulation{
			string compute[3], geometry;
			double dx, dy, dz, dv;
			double nSets, nEquilSets, nStepsPerSet, printEvery;
			double displaceProb, exchangeProb, volumeProb;
			int step;
			bool PBC[3];
		} sim;
		struct Fluid{
			string name, vdwPot;
			double temp, dens, molarMass, extPress, boxWidth, mu;
			double rcut, sigma, epsilon;
			double systemEnergy;
			int nParts;
		} fluid;
		struct Particle{
			double x, y, z;
		} part[MAXPART];
/* **************************************************************************** */
};
/* **************************************************************************** */

#endif /* _MC_H_ */

