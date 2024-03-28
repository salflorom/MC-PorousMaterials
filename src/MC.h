//Author: Santiago A. Flores Roman

#ifndef _MC_H_
#define _MC_H_

#include <string> // string
#include <cmath> //M_PI
#include <random> // random_device, mt19937_64, uniform_real_distribution
#include "Eigen/Dense"

#define NBINS 100
#define MAXPART 1000000 //Max. num. of particles. 0th part. is to save old config.
#define MAXSPECIES 10 //Max. num. of species.
#define MAXBOX 10 //Max. num. of boxes. 0th box is to save old comfiguration.
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
			Eigen::Matrix<int,MAXBOX,1> acceptance, rejection, nDisplacements;
			Eigen::Matrix<int,MAXBOX,1> acceptanceVol, rejectionVol;
			int nVolChanges;
			int acceptSwap, rejectSwap, nSwaps;
			Eigen::Matrix<int,MAXBOX,MAXSPECIES> widomInsertions, widomDeletions;
			Eigen::Matrix<double,MAXBOX,MAXSPECIES> barC, widomInsert, widomDelete;
			Eigen::Matrix<double,Eigen::Dynamic,1> rdf = Eigen::Matrix<double,Eigen::Dynamic,1>(NBINS+1);
			double binWidth;
		} stats;
		struct Simulation{
			bool printTrajectory, continueAfterCrash;
			double nSets, nEquilSets, nStepsPerSet, printEvery;
			Eigen::Matrix<double,MAXBOX,1> dr, dv;
			Eigen::Vector2i rdf;
			int nDispAttempts, nSwapAttempts, nVolAttempts, cycle;
			string projName;
		} sim;
		struct Particle{
			double x, y, z;
			double energy, manyBodyE, pairPotE, boxE; // Particle energy.
			bool operator!=(Particle bPart){
				if ((this->x != bPart.x) || (this->y != bPart.y) || (this->z != bPart.z)){
					return true;
				}else return false;
			}
		};
		struct Fluid{
			string name;
			Eigen::Matrix<string,MAXSPECIES,1> vdwPot;
			int nParts;
			double molarMass, mu, muEx;
			Eigen::Matrix<double,MAXSPECIES,1> sigma, epsilon, rcut;
			Eigen::Matrix<Particle,Eigen::Dynamic,1> particle = Eigen::Matrix<Particle,Eigen::Dynamic,1>(MAXPART);
		};
		Eigen::Matrix<Fluid,Eigen::Dynamic,1> fluid = Eigen::Matrix<Fluid,Eigen::Dynamic,1>(MAXSPECIES); //For single species properties.
		struct Box{
			string name, geometry, vdwPot; //pot: solid-fluid potential.
			int nParts, fix;
			Eigen::Vector3i PBC;
			Eigen::Vector3d width;
			double oldEnergy, solidDens, volume, deltaLayers, maxRcut;
			double energy, manyBodyE, pairPotE, boxE; // Energy of the box.
			int nLayersPerWall;
			Eigen::Matrix<Fluid,Eigen::Dynamic,1> fluid = Eigen::Matrix<Fluid,Eigen::Dynamic,1>(MAXSPECIES); //For single species properties.
		};
		Eigen::Matrix<Box,Eigen::Dynamic,1> box = Eigen::Matrix<Box,Eigen::Dynamic,1>(MAXBOX); //For single species properties.
		struct ThermodynamicSystem{
			double temp, volume, press;
			int nBoxes, nSpecies, nParts;
		} thermoSys;
/* **************************************************************************** */
	public:
/* **************************************************************************** */
		MC(void);
		void OutputDirectory(void);
		size_t ReadInputFile(string);
		void ResetStats(void);
		void PrintParams(void);
		void InsertParticle(int, int, int);
		void InitialConfig(void);
		void MinimizeEnergy(size_t);
		void AdjustMCMoves(size_t);
		void MoveParticle(void);
		void ChangeVolume(void);
		void RescaleCenterOfMass(Box&, Box&, int);
		void SwapParticle(void);
		void PBC(int, Particle&);
		void ComputeRDF(void);
		void PrintRDF(void);
		void Metropolis(int);
		void PrintTrajectory(size_t);
		void PrintLog(size_t);
		void PrintStats(size_t);
		void ComputeWidom(void);
		void ComputeChemicalPotential(void);
		void EnergyOfParticle(int, int, int);
		void BoxEnergy(int);
		void CorrectEnergy(void);

		void ReadTrajectory(void);
		size_t ReadLogFile(void);
		size_t GetNSets(void){return int(sim.nSets);}
		size_t GetNEquilSets(void){return int(sim.nEquilSets);}
		size_t GetNStepsPerSet(void){return int(sim.nStepsPerSet);}
		size_t GetPrintEvery(void){return int(sim.printEvery);}
		Eigen::Vector3i GetMCMoves(void);

		double Random(void){return dis(engine);} //random num. in the interval [0,1).
		double BarWeight(double, double);
		double ComputeVolume(Box&);
		double ComputeBoxWidth(Box&, double);

		double NeighDistance(int, Particle, Particle);
		// >>> Fluid-Fluid potentials >>> //
		// Hard sphere potential
		double HardSphere_Pot(int, int, int, int);
		// Lennard-Jones 12-6 potential
		double LJ126_Pot(int, int, int, int);
		// EAM Na potential
		Eigen::Vector2d EAMNa_Pot(int, int, int, int);
		double EAMNa_EmbPot(double);
		double EAMNa_eDens(double);
		double EAMNa_PairPot(double);
		// EAM K potential
		Eigen::Vector2d EAMK_Pot(int, int, int, int);
		double EAMK_EmbPot(double);
		double EAMK_eDens(double);
		double EAMK_PairPot(double);
		// EAM Rb potential
		Eigen::Vector2d EAMRb_Pot(int, int, int, int);
		double EAMRb_EmbPot(double);
		double EAMRb_eDens(double);
		double EAMRb_PairPot(double);
		// EAM Ga potential
		Eigen::Vector2d EAMGa_Pot(int, int, int, int);
		double StepUnit(double, double, double);
		double EAMGa_EmbPot(double);
		double EAMGa_eDens(double);
		double EAMGa_PairPot(double);
		// <<< Fluid-Fluid potentials <<< //

		// >>> Solid-Fluid potentials >>> //
		double SlitLJ_Pot(int, int, int);
		double CylindricalLJ10_4(int, int, int);
		double CylindricalSteele10_4_3(int, int, int);
		double SphericalLJ_Pot(int, int, int);
		// <<< Solid-Fluid potentials >>> //
/* **************************************************************************** */
};
/* **************************************************************************** */

#endif /* _MC_H_ */

