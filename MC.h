//Author: Santiago A. Flores Roman

#ifndef _MC_H_
#define _MC_H_

#include <string>
#include <cstdlib>
#include <cmath>
#include <vector>
#include <random>

#define NBINS 500
#define MAXPART 99999 //Max. num. of particles.
#define pi M_PI
#define kb 1.380649e-23 // J/K
#define na 6.02214076e23 // mol^-1
#define planck 6.62607015e-34 // J*s

using namespace std;

/* **************************************************************************** */
class MC {
/* **************************************************************************** */
	public:
/* **************************************************************************** */
		MC(void);
		string OutputDirectory(bool createDir);
		void ReadInputFile(string inFileName);
		void ResetStats(void);
		void PrintParams(void);
		void InsertParticle(int);
		void InitialConfig(void);
		void MinimizeEnergy(void);
		void AdjustMCMoves(void);
		void MoveParticle(void);
		void ChangeVolume(void);
		void RescaleCenterOfMass(double initBoxWidth, double finalBoxWidth);
		void ExchangeParticle(void);
		void PBC(int index);
		void ComputeRDF(void);
		void PrintRDF(void);
		void Metropolis(int index);
		void CreateEXYZ(int step);
		void CreateLogFile(int step);
		void PrintStats(int step);
		void ComputeWidom(void);
		void ComputeChemicalPotential(void);
		void EnergyOfParticle(int index);
		double Random(void){return dis(engine);} //random num. in the interval [0,1).
		double ComputeVolume(void);
		double ComputeBoxWidth(double volume);
		double* SystemEnergy(void);
		double NeighDistance(int i, int j);

		// Fluid-Fluid potentials //
		double LJ_Energy(int i, int j);
		// EAM Ga potential vvvvv //
		double* EAMGA_Energy(int index);
		double StepUnit(double radius, double leftLim, double rightLim);
		double EmbPot(double rho);
		double eDens(double radius);
		double PairPot(double radius);
		// EAM Ga potential ^^^^^ //
		// Fluid-Fluid potentials //

		// Solid-Fluid potentials //
		double SlitLJ(int index);
		double CylindricalLJ(int index);
		double SphericalLJ(int index);
		// Solid-Fluid potentials //

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
			string rdf, compute[3];
			double dr, dv;
			double nSets, nEquilSets, nStepsPerSet, printEvery;
			double displaceProb, exchangeProb, volumeProb;
			double systemEnergy, deltaParticleEnergy;
			int step;
		} sim;
		struct FrameWork{
			string name, sfPot, geometry;
			bool PBC[3];
			double boxWidth[3];
			double sfEnergy, deltasfEnergy;
			double sfDensity, sfEpsilon, sfSigma, deltaLayers;
			int nLayersPerWall;
		} pore;
		struct Fluid{
			string name, vdwPot;
			double temp, molarMass, dens, extPress, mu;
			double rcut, sigma, epsilon;
			double ffEnergy, deltaffManybody, deltaffPairPot;
			int nParts;
		} fluid;
		struct Particle{
			double x, y, z;
		} part[MAXPART];
/* **************************************************************************** */
};
/* **************************************************************************** */

#endif /* _MC_H_ */

