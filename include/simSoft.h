//
// Author: Francesco Arceri
// Date: March 28 2025
//
// HEADER FILE FOR SIMSOFT CLASS

#ifndef _H
#define SIMSOFT_H

#include "defs.h"
#include <vector>
#include <tuple>
#include <string>
#include <memory>
#include <iomanip>

using namespace std;
using std::vector;
using std::string;
using std::tuple;

struct simControlStruct {
  enum class particleEnum {passive, active, vicsek} particleType;
  enum class boundaryEnum {pbc, fixed} boundaryType;
  enum class neighborEnum {neighbor, allToAll} neighborType;
  enum class potentialEnum {harmonic, lennardJones, WCA, doubleLJ, LJMinusPlus, LJWCA} potentialType;
};

// pointer-to-member function call macro
#define CALL_MEMBER_FN(object, ptrToMember) ((object).*(ptrToMember))

class simSoft;
class FIRE;
class SimInterface;

class simSoft
{
public:

  // constructor and deconstructor
  simSoft(long nParticles, long dim);
  ~simSoft();

  // Simulator
  FIRE * fire_;
  SimInterface * sim_;

  simControlStruct simControl;

  // sp packing constants
  long nDim;
  long numParticles;

  std::vector<double> boxSize;

  // time step
  double dt;
  // energy scale
  double ec;
  // self-propulsion parameters
  double driving, taup;
  // Vicsek velocity interaction parameters
  double Jvicsek, Rvicsek;
  // Lennard-Jones constants
  double LJcutoff, LJecut, LJfshift;
  double LJecutPlus, LJfshiftPlus;
  // double Lennard-Jones constants
  double eAA, eAB, eBB;
  long numA;
  // neighbor update variables
  double cutoff, cutDistance, rcut;
  long updateCount;
  bool shift;

  // dynamical particle variables
  std::vector<double> rad;
  std::vector<double> pos;
  std::vector<double> vel;
  std::vector<double> squaredVel;
  std::vector<double> force;
  std::vector<double> energy;
  std::vector<double> angle;
  std::vector<double> omega;
  std::vector<double> alpha;
  std::vector<double> randAngle;
  std::vector<double> randomAngle;

  // correlation variables
  std::vector<double> initPos;
  std::vector<double> lastPos;
  std::vector<double> vicsekLastPos;
	std::vector<double> delta;
  std::vector<double> disp;

	// neighbor list
  std::vector<int> dispFlag;
  std::vector<long> neighborList;
  std::vector<long> maxNeighborList;
  long maxNeighbors;
	long neighborListSize;
  long neighborLimit;

  // Vicsek interaction list
  std::vector<int> vicsekFlag;
  std::vector<long> vicsekNeighborList;
  std::vector<long> vicsekMaxNeighborList;
  long vicsekMaxNeighbors;
	long vicsekNeighborListSize;
  long vicsekNeighborLimit;

  void initVariables(long numParticles_);

  void initDeltaVariables(long numParticles_);

  void initNeighbors(long numParticles_);

  void initVicsekNeighbors(long numParticles_);

  //setters and getters
  void setParticleType(simControlStruct::particleEnum particleType_);

  void setBoundaryType(simControlStruct::boundaryEnum boundaryType_);
	simControlStruct::boundaryEnum getBoundaryType();

  void setNeighborType(simControlStruct::neighborEnum neighborType_);
	simControlStruct::neighborEnum getNeighborType();

  void setPotentialType(simControlStruct::potentialEnum potentialType_);
	simControlStruct::potentialEnum getPotentialType();

  void setNDim(long nDim_);
  long getNDim();

  void setNumParticles(long numParticles_);
	long getNumParticles();

  long getTypeNumParticles();

  void setBoxSize(std::vector<double> &boxSize_);
  std::vector<double> getBoxSize();

  void setRadii(std::vector<double> &rad_);
  std::vector<double> getRadii();

  double getMeanParticleSize();

  void checkPBC();

  void setPositions(std::vector<double> &pos_);
  std::vector<double> getPositions();

  void setPBCPositions(std::vector<double> &pos_);
  std::vector<double> getPBCPositions();

  void resetLastPositions();

  void resetVicsekLastPositions();

  void setInitialPositions();

  std::vector<double> getLastPositions();

  void setVelocities(std::vector<double> &vel_);
  std::vector<double> getVelocities();

  void setForces(std::vector<double> &force_);
  std::vector<double> getForces();

  std::vector<double> getEnergies();

  void setAngles(std::vector<double> &angle_);
  std::vector<double> getAngles();

  double getPackingFraction();

  // initialization functions
  void setPolyRandomParticles(double phi0, double polyDispersity, double lx, double ly, double lz);

  void setMonoRandomParticles(double phi0, double lx, double ly, double lz);

  void setBiRandomParticles(double phi0, double lx, double ly, double lz);

  void scaleParticles(double scale);

  void scaleParticlePacking();

  void scaleVelocities(double scale);

  void initializeParticleAngles();

  // force and energy
  void setEnergyCostant(double ec_);
  double getEnergyCostant();

  double setTimeStep(double dt_);

  void setSelfPropulsionParams(double driving_, double taup_);
  void getSelfPropulsionParams(double &driving_, double &taup_);

  void setVicsekParams(double driving_, double taup_, double Jvicsek_, double Rvicsek_);
  void getVicsekParams(double &driving_, double &taup_, double &Jvicsek_, double &Rvicsek_);

  void setLJcutoff(double LJcutoff_);

  void setDoubleLJconstants(double LJcutoff_, double eAA_, double eAB_, double eBB_, long numA_);

  void setLJMinusPlusParams(double LJcutoff_, long numA_);

  void setLJWCAparams(double LJcutoff_, long numA_);

  void addSelfPropulsion();

  void addVicsekAlignment();

  void calcVicsekAlignment();

  void calcInteraction();

  void addFixedWallInteraction();

  void calcForceEnergy();

  double getMaxUnbalancedForce();

  double getPotentialEnergy();

  double getKineticEnergy();

  double getTemperature();

  double getEnergy();

  void adjustKineticEnergy(double prevEtot);

  // neighbor update
  double setDisplacementCutoff(double cutoff_);

  void resetUpdateCount();

  long getUpdateCount();

  void removeCOMDrift();

  void checkParticleDisplacement();

  void checkParticleNeighbors();

  void checkVicsekNeighbors();

  std::vector<long> getNeighbors();

  std::vector<long> getVicsekNeighbors();

  void calcNeighbors(double cutDistance);

  void calcNeighborList(double cutDistance);

  void calcVicsekNeighborList();

  // minimizer functions
  void initFIRE(std::vector<double> &FIREparams, long minStep, long maxStep, long numDOF);

  void setParticleMassFIRE();

  void setTimeStepFIRE(double timeStep);

  void FIRELoop();

  // Langevin integrator
  void initSoftLangevin(double Temp, double gamma, bool readState);

  void softLangevinLoop(bool conserve = false);

  // NVE integrator
  void initSoftNVE(double Temp, bool readState);

  void softNVELoop();

  void rescaleVelocities(double Temp);
  
  // Nose-Hoover integrator
  void getNoseHooverParams(double &mass, double &damping);

  void initSoftNoseHoover(double Temp, double mass, double gamma, bool readState);

  void softNoseHooverLoop();

  // simple test functions
  void setTwoParticleTestPacking(double sigma0, double sigma1, double lx, double ly, double y0, double y1, double vel1);

  void setThreeParticleTestPacking(double sigma01, double sigma2, double lx, double ly, double y01, double y2, double vel2);

  void firstUpdate(double timeStep);

  void secondUpdate(double timeStep);

  void testInteraction(double timeStep);

};

#endif /* SP2D_H */
