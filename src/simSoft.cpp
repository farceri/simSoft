//
// Author: Francesco Arceri
// Date: March 28 2025
//
// FUNCTION DECLARATIONS

#include "../include/simSoft.h"
#include "../include/defs.h"
#include "../include/Integrator.h"
#include "../include/FIRE.h"
#include <iostream>
#include <string>
#include <vector>
#include <cmath>
#include <random>
#include <time.h>
#include <cuda.h>

using namespace std;
using std::cout;
using std::endl;

//************************** sp object definition ***************************//
simSoft::simSoft(long nParticles, long dim) {
  // default values
  srand48(time(0));
  nDim = dim;
  numParticles = nParticles;
  numA = numParticles;
  setNDim(nDim);
  setNumParticles(numParticles);
	simControl.particleType = simControlStruct::particleEnum::passive;
	simControl.boundaryType = simControlStruct::boundaryEnum::pbc;
	simControl.neighborType = simControlStruct::neighborEnum::neighbor;
	simControl.potentialType = simControlStruct::potentialEnum::harmonic;
  // default parameters
  dt = 1e-04;
	ec = 1;
  cutDistance = 1;
  updateCount = 0;
  shift = false;
  d_boxSize.resize(nDim);
  std::fill(boxSize.begin(), boxSize.end(), double(1));
  // initialize vectors
  initVariables(numParticles);
  initDeltaVariables(numParticles);
  initNeighbors(numParticles);
}

simSoft::~simSoft() {}

void simSoft::initVariables(long numParticles_) {
  rad.resize(numParticles_);
  pos.resize(numParticles_ * nDim);
  vel.resize(numParticles_ * nDim);
  force.resize(numParticles_ * nDim);
  energy.resize(numParticles_);
  squaredVel.resize(numParticles_ * nDim);
  std::fill(rad.begin(), rad.end(), double(0));
  std::fill(pos.begin(), pos.end(), double(0));
  std::fill(vel.begin(), vel.end(), double(0));
  std::fill(force.begin(), force.end(), double(0));
  std::fill(energy.begin(), energy.end(), double(0));
  std::fill(squaredVel.begin(), squaredVel.end(), double(0));
}

void simSoft::initDeltaVariables(long numParticles_) {
  initPos.resize(numParticles_ * nDim);
  lastPos.resize(numParticles * nDim);
  delta.resize(numParticles_ * nDim);
  disp.resize(numParticles_);
  std::fill(initPos.begin(), initPos.end(), double(0));
  std::fill(lastPos.begin(), lastPos.end(), double(0));
  std::fill(delta.begin(), delta.end(), double(0));
  std::fill(disp.begin(), disp.end(), double(0));
}

void simSoft::initNeighbors(long numParticles_) {
  neighborListSize = 0;
  maxNeighbors = 0;
  neighborList.resize(numParticles_);
  maxNeighborList.resize(numParticles_);
  std::fill(neighborList.begin(), neighborList.end(), -1L);
  std::fill(maxNeighborList.begin(), maxNeighborList.end(), maxNeighbors);
  dispFlag.resize(numParticles_);
  std::fill(dispFlag.begin(), dispFlag.end(), int(0));
}

void simSoft::initVicsekNeighbors(long numParticles_) {
  vicsekNeighborListSize = 0;
  vicsekMaxNeighbors = 0;
  vicsekNeighborList.resize(numParticles_);
  vicsekMaxNeighborList.resize(numParticles_);
  std::fill(vicsekNeighborList.begin(), vicsekNeighborList.end(), -1L);
  std::fill(vicsekMaxNeighborList.begin(), vicsekMaxNeighborList.end(), vicsekMaxNeighbors);
  vicsekFlag.resize(numParticles_);
  std::fill(vicsekFlag.begin(), vicsekFlag.end(), int(0));
}

//**************************** setters and getters ***************************//
void simSoft::setParticleType(simControlStruct::particleEnum particleType_) {
	simControl.particleType = particleType_;
  if(simControl.particleType == simControlStruct::particleEnum::passive) {
    cout << "simSoft::setParticleType: particleType: passive" << endl;
  } else if(simControl.particleType == simControlStruct::particleEnum::active) {
    if(nDim == 2) {
      angle.resize(numParticles);
      randAngle.resize(numParticles);
    } else if(nDim == 3) {
      angle.resize(numParticles * nDim);
      randAngle.resize(numParticles * nDim);
    }
    std::fill(angle.begin(), angle.end(), double(0));
    std::fill(randAngle.begin(), randAngle.end(), double(0));
    cout << "simSoft::setParticleType: particleType: active" << endl;
  } else if(simControl.particleType == simControlStruct::particleEnum::vicsek) {
    randAngle.resize(numParticles);
    angle.resize(numParticles);
    alpha.resize(numParticles);
    std::fill(randAngle.begin(), randAngle.end(), double(0));
    std::fill(angle.begin(), angle.end(), double(0));
    std::fill(alpha.begin(), alpha.end(), double(0));
    initVicsekNeighbors(numParticles);
    d_vicsekLastPos.resize(numParticles * nDim);
    std::fill(d_vicsekLastPos.begin(), d_vicsekLastPos.end(), double(0));
    cout << "simSoft::setParticleType: particleType: vicsek" << endl;
  } else {
    cout << "simSoft::setParticleType: please specify valid particleType: passive, active or vicsek" << endl;
  }
}

void simSoft::setBoundaryType(simControlStruct::boundaryEnum boundaryType_) {
	simControl.boundaryType = boundaryType_;
  if(simControl.boundaryType == simControlStruct::boundaryEnum::pbc) {
    cout << "simSoft::setBoundaryType: boundaryType: pbc" << endl;
  } else if(simControl.boundaryType == simControlStruct::boundaryEnum::fixed) {
    cout << "simSoft::setBoundaryType: boundaryType: fixed" << endl;
  } else {
    cout << "simSoft::setBoundaryType: please specify valid boundaryType: pbc, leesEdwards, fixed, reflect, reflectNoise, rigid, mobile and plastic" << endl;
  }
}

simControlStruct::boundaryEnum simSoft::getBoundaryType() {
	return simControl.boundaryType;
}

void simSoft::setNeighborType(simControlStruct::neighborEnum neighborType_) {
	simControl.neighborType = neighborType_;
  if(simControl.neighborType == simControlStruct::neighborEnum::neighbor) {
    cout << "simSoft::setNeighborType: neighborType: neighbor" << endl;
  } else if(simControl.neighborType == simControlStruct::neighborEnum::allToAll) {
    cout << "simSoft::setNeighborType: neighborType: allToAll" << endl;
  } else {
    cout << "simSoft::setNeighborType: please specify valid neighborType: neighbor or allToAll" << endl;
  }
}

simControlStruct::neighborEnum simSoft::getNeighborType() {
	return simControl.neighborType;
}

void simSoft::setPotentialType(simControlStruct::potentialEnum potentialType_) {
	simControl.potentialType = potentialType_;
  if(simControl.potentialType == simControlStruct::potentialEnum::harmonic) {
    cout << "simSoft::setPotentialType: potentialType: harmonic" << endl;
  } else if(simControl.potentialType == simControlStruct::potentialEnum::lennardJones) {
    cout << "simSoft::setPotentialType: potentialType: lennardJones" << endl;
  } else if(simControl.potentialType == simControlStruct::potentialEnum::Mie) {
    cout << "simSoft::setPotentialType: potentialType: Mie" << endl;
  } else if(simControl.potentialType == simControlStruct::potentialEnum::WCA) {
    cout << "simSoft::setPotentialType: potentialType: WCA" << endl;
  } else if(simControl.potentialType == simControlStruct::potentialEnum::doubleLJ) {
    cout << "simSoft::setPotentialType: potentialType: doubleLJ" << endl;
  } else if(simControl.potentialType == simControlStruct::potentialEnum::LJMinusPlus) {
    cout << "simSoft::setPotentialType: potentialType: LJMinusPlus" << endl;
  } else if(simControl.potentialType == simControlStruct::potentialEnum::LJWCA) {
    cout << "simSoft::setPotentialType: potentialType: LJWCA" << endl;
  } else {
    cout << "simSoft::setPotentialType: please specify valid potentialType: none, harmonic, lennardJones, WCA, adhesive, doubleLJ, LJMinusPlus and LJWCA" << endl;
  }
}

simControlStruct::potentialEnum simSoft::getPotentialType() {
	return simControl.potentialType;
}

void simSoft::setNDim(long nDim_) {
  nDim = nDim_;
}

long simSoft::getNDim() {
	return nDim;
}

void simSoft::setNumParticles(long numParticles_) {
  numParticles = numParticles_;
}

long simSoft::getNumParticles() {
	return numParticles;
}

long simSoft::getTypeNumParticles() {
	return numA;
}

//TODO: error messages for all the vector getters and setters
void simSoft::setBoxSize(std::vector<double> &boxSize_) {
  if(boxSize_.size() == ulong(nDim)) {
    boxSize = boxSize_;
  } else {
    cout << "simSoft::setBoxSize: size of boxSize_ does not match nDim" << endl;
  }
}

std::vector<double> simSoft::getBoxSize() {
  return boxSize;
}

void simSoft::setRadii(std::vector<double> &rad_) {
  if(rad_.size() == ulong(numParticles)) {
    rad = rad_;
  } else {
    cout << "simSoft::setRadii: size of rad_ does not match numParticles" << endl;
  }
}

std::vector<double> simSoft::getRadii() {
  return rad;
}

double simSoft::getMeanParticleSize() {
  return 2 * std::reduce(rad.begin(), rad.end(), double(0), std::plus<double>()) / numParticles;
}

void simSoft::setPositions(std::vector<double> &pos_) {
  if(pos_.size() == ulong(numParticles * nDim)) {
    pos = pos_;
  } else {
    cout << "simSoft::setPositions: size of pos_ does not match numParticles * nDim" << endl;
  }
}

std::vector<double> simSoft::getPositions() {
  return pos;
}

void simSoft::checkPBC() {
  std::vector<double> posPBC(pos.size());
  const double *pPos = pos.data();
  double *pPosPBC = posPBC.data();
  double *pBoxSize = boxSize.data();

  long f_nDim(nDim);
  auto checkPBC = [&](long pId) {
    for (long dim = 0; dim < f_nDim; dim++) {
			pPosPBC[pId * f_nDim + dim] = pPos[pId * f_nDim + dim] - floor(pPos[pId * f_nDim + dim] / pBoxSize[dim]) * pBoxSize[dim];
		}
  };

	std::for_each(r, r + vel.size(), checkPBC);
  pos = posPBC;
}

void simSoft::setPBCPositions(std::vector<double> &pos_) {
  if(pos_.size() == ulong(numParticles * nDim)) {
    pos = pos_;
  } else {
    cout << "simSoft::setPBCPositions: size of pos_ does not match numParticles * nDim" << endl;
  }
  checkPBC();
}


std::vector<double> simSoft::getPBCPositions() {
  checkPBC();
  return pos;
}

void simSoft::resetLastPositions() {
  lastPos = pos;
}

void simSoft::resetVicsekLastPositions() {
  vicsekLastPos = pos;
}

void simSoft::setInitialPositions() {
  initPos = pos;
}

std::vector<double> simSoft::getLastPositions() {
  return lastPos;
}

void simSoft::setVelocities(std::vector<double> &vel_) {
  if(vel_.size() == ulong(numParticles * nDim)) {
    vel = vel_;
  } else {
    cout << "simSoft::setVelocities: size of vel_ does not match numParticles * nDim" << endl;
  }
}

std::vector<double> simSoft::getParticleVelocities() {
  return vel;
}

void simSoft::setForces(std::vector<double> &force_) {
  if(force_.size() == ulong(numParticles * nDim)) {
    force = force_;
  } else {
    cout << "simSoft::setForces: size of force_ does not match numParticles * nDim" << endl;
  }
}

std::vector<double> simSoft::getParticleForces() {
  return force;
}

std::vector<double> simSoft::getParticleEnergies() {
  return energy;
}

void simSoft::setAngles(std::vector<double> &angle_) {
  if(nDim == 2) {
    if(angle_.size() == ulong(numParticles)) {
      angle = angle_;
    } else {
      cout << "simSoft::setAngles: size of force_ does not match numParticles * nDim" << endl;
    }
  } else if(nDim == 3) {
    if(angle_.size() == ulong(numParticles * nDim)) {
      angle = angle_;
    } else {
      cout << "simSoft::setAngles: size of force_ does not match numParticles * nDim" << endl;
    }
  } else {
    cout << "simSoft::setAngles: only dimensions 2 and 3 are allowed!" << endl;
  }
}

std::vector<double> simSoft::getAngles() {
  return angle;
}

double simSoft::getPackingFraction() {
  if(nDim == 2) {
    std::vector<double> radSquared(numParticles);
    std::transform(rad.begin(), rad.end(), radSquared.begin(), square());
    return std::reduce(radSquared.begin(), radSquared.end(), double(0), std::plus<double>()) * PI / (boxSize[0] * boxSize[1]);
  } else if(nDim == 3) {
    std::vector<double> radCubed(numParticles);
    std::transform(rad.begin(), rad.end(), radCubed.begin(), cube());
    return std::reduce(radCubed.begin(), radCubed.end(), double(0), std::plus<double>()) * 4 * PI / (3 * boxSize[0] * boxSize[1] * boxSize[2]);
  } else {
    cout << "simSoft::getPackingFraction: only dimensions 2 and 3 are allowed!" << endl;
    return 0;
  }
}

//************************ initialization functions **************************//
vid simSoft::setPolyRandomParticles(double phi0, double polyDispersity, double lx, double ly, double lz) {
  double r1, r2, randNum, mean = 0, sigma, scale;
  sigma = sqrt(log(polyDispersity*polyDispersity + 1.));
  // generate polydisperse particle size
  for (long particleId = 0; particleId < numParticles; particleId++) {
    r1 = drand48();
    r2 = drand48();
    randNum = sqrt(-2. * log(r1)) * cos(2. * PI * r2);
    rad[particleId] = 0.5 * exp(mean + randNum * sigma);
  }
  boxSize[0] = lx;
  boxSize[1] = ly;
  if(nDim == 2) {
    scale = sqrt(getParticlePhi() / phi0);
  } else if(nDim == 3) {
    boxSize[2] = lz;
    scale = cbrt(getParticlePhi() / phi0);
    boxSize[2] = lz * scale;
  } else {
    cout << "simSoft::setScaledPolyRandomSoftParticles: only dimesions 2 and 3 are allowed!" << endl;
  }
  boxSize[0] = lx * scale;
  boxSize[1] = ly * scale;
  // extract random positions
  for (long particleId = 0; particleId < numParticles; particleId++) {
    for(long dim = 0; dim < nDim; dim++) {
      pos[particleId * nDim + dim] = boxSize[dim] * drand48();
    }
  }
}

void simSoft::setMonoRandomParticles(double phi0, double lx, double ly, double lz) {
  double scale;
  // generate polydisperse particle size
  std::fill(d_particleRad.begin(), d_particleRad.end(), 0.5);
  boxSize[0] = lx;
  boxSize[1] = ly;
  if(nDim == 2) {
    scale = sqrt(getParticlePhi() / phi0);
  } else if(nDim == 3) {
    boxSize[2] = lz;
    scale = cbrt(getParticlePhi() / phi0);
    boxSize[2] = lz * scale;
  } else {
    cout << "simSoft::setScaledPolyRandomSoftParticles: only dimesions 2 and 3 are allowed!" << endl;
  }
  boxSize[0] = lx * scale;
  boxSize[1] = ly * scale;
  // extract random positions
  for (long particleId = 0; particleId < numParticles; particleId++) {
    for(long dim = 0; dim < nDim; dim++) {
      pos[particleId * nDim + dim] = boxSize[dim] * drand48();
    }
  }
}

void simSoft::setBiRandomParticles(double phi0, double lx, double ly, double lz) {
  double scale;
  long halfNum = int(numParticles / 2);
  // generate polydisperse particle size
  std::fill(d_particleRad.begin(), d_particleRad.begin() + halfNum, 0.5);
  std::fill(d_particleRad.begin() + halfNum, d_particleRad.end(), 0.7);
  boxSize[0] = lx;
  boxSize[1] = ly;
  if(nDim == 2) {
    scale = sqrt(getParticlePhi() / phi0);
  } else if(nDim == 3) {
    boxSize[2] = lz;
    scale = cbrt(getParticlePhi() / phi0);
    boxSize[2] = lz * scale;
  } else {
    cout << "simSoft::setScaledBiRandomSoftParticles: only dimesions 2 and 3 are allowed!" << endl;
  }
  boxSize[0] = lx * scale;
  boxSize[1] = ly * scale;
  // extract random positions
  for (long particleId = 0; particleId < numParticles; particleId++) {
    for(long dim = 0; dim < nDim; dim++) {
      pos[particleId * nDim + dim] = boxSize[dim] * drand48();
    }
  }
}

void simSoft::scaleParticles(double scale) {
  std::transform(rad.begin(), rad.end(), rad.begin(), [](double x) { return x * scale; });
}

void simSoft::scaleParticlePacking() {
  double sigma = getMeanParticleSize();
  std::transform(rad.begin(), rad.end(), rad.begin(), [](double x) { return x / scale; });
  std::transform(pos.begin(), pos.end(), pos.begin(), [](double x) { return x / scale; });
  for (long dim = 0; dim < nDim; dim++) {
    boxSize[dim] /= sigma;
  }
}

void simSoft::scaleVelocities(double scale) {
  std::transform(vel.begin(), vel.end(), vel.begin(), [](double x) { return x * scale; });
}

// compute particle angles from velocity
void simSoft::initializeParticleAngles() {
  long f_nDim(nDim);
  randNum generateRandomNumbers(-PI, PI);  // Create a randNum object with the desired range
  std::for_each(angle.begin(), angle.end(), [&](double& val) {val = generateRandomNumbers(0);});
  if(nDim == 3) {
      double* pAngle = angle.data();
      const double* pVel = vel.data();

      auto compute3DAngle = [&](long particleId) {
      auto theta = acos(pVel[particleId * f_nDim + 2]);
      auto phi = atan(pVel[particleId * f_nDim + 1] / pVel[particleId * f_nDim]);
      pAngle[particleId * f_nDim] = cos(theta) * cos(phi);
      pAngle[particleId * f_nDim + 1] = sin(theta) * cos(phi);
      pAngle[particleId * f_nDim + 2] = sin(phi);
    };

    std::for_each(r, r + numParticles, compute3DAngle);
  }
}

//*************************** force and energy *******************************//
void simSoft::setEnergyCostant(double ec_) {
  ec = ec_;
}

double simSoft::getEnergyCostant() {
  if(simControl.potentialType == simControlStruct::potentialEnum::doubleLJ) {
    return sqrt(eAA * eBB);
  } else {
    return ec;
  }
}

double simSoft::setTimeStep(double dt_) {
  dt = dt_;
  return dt;
}

void simSoft::setSelfPropulsionParams(double driving_, double taup_) {
  driving = driving_;
  taup = taup_;
  cout << "simSoft::setSelfPropulsionParams:: driving: " << driving << " taup: " << taup << endl;
}

void simSoft::getSelfPropulsionParams(double &driving_, double &taup_) {
  driving_ = driving;
  taup_ = taup;
  //cout << "simSoft::getSelfPropulsionParams:: driving: " << driving_ << " taup: " << taup_ << endl;
}

void simSoft::setVicsekParams(double driving_, double taup_, double Jvicsek_, double Rvicsek_) {
  driving = driving_;
  taup = taup_;
  Jvicsek = Jvicsek_;
  Rvicsek = Rvicsek_;
  cout << "simSoft::setVicsekParams:: driving: " << driving << " interactin strength: " << Jvicsek << " and radius: " << Rvicsek << endl;
}

void simSoft::getVicsekParams(double &driving_, double &taup_, double &Jvicsek_, double &Rvicsek_) {
  driving_ = driving;
  taup_ = taup;
  Jvicsek_ = Jvicsek;
  Rvicsek_ = Rvicsek;
  //cout << "simSoft::getVicsekParams:: driving: " << driving_ << " interactin strength: " << Jvicsek << " and radius: " << Rvicsek << endl;
}

void simSoft::setLJcutoff(double LJcutoff_) {
  LJcutoff = LJcutoff_;
  double ratio6 = 1 / pow(LJcutoff, 6);
  LJecut = 4 * ec * (ratio6 * ratio6 - ratio6);
	LJfshift = 24 * ec * (2 * ratio6 - 1) * ratio6 / LJcutoff;
  cout << "simSoft::setLJcutoff::LJcutoff: " << LJcutoff << " energy shift: " << LJecut << " LJfshift: " << LJfshift << endl;
}

void simSoft::setDoubleLJconstants(double LJcutoff_, double eAA_, double eAB_, double eBB_, long numA_) {
  LJcutoff = LJcutoff_;
  eAA = eAA_;
  eAB = eAB_;
  eBB = eBB_;
  double ratio6 = 1 / pow(LJcutoff, 6);
  LJecut = 4 * (ratio6 * ratio6 - ratio6);
	LJfshift = 24 * (2 * ratio6 - 1) * ratio6 / LJcutoff;
  numA = numA_;
  cout << "simSoft::setDoubleLJconstants::eAA: " << eAA << " eAB: " << eAB << " eBB: " << eBB;
  cout << " LJcutoff: " << LJcutoff << " LJecut: " << LJecut << " LJfshift: " << LJfshift;
  cout << " numA: " << numA << endl;
}

void simSoft::setLJMinusPlusParams(double LJcutoff_, long numA_) {
  LJcutoff = LJcutoff_;
  double ratio6 = 1 / pow(LJcutoff, 6);
  LJecut = 4 * ec * (ratio6 * ratio6 - ratio6);
	LJfshift = 24 * ec * (2 * ratio6 - 1) * ratio6 / LJcutoff;
  // repulsive Lennard-Jones
  LJecutPlus = 4 * ec * (ratio6 * ratio6 + ratio6);
	LJfshiftPlus = 24 * ec * (2 * ratio6 + 1) * ratio6 / LJcutoff;
  cout << "simSoft::setLJMinusPlusParams::LJcutoff: " << LJcutoff << " energy shift: " << LJecut << " LJfshift: " << LJfshift << endl;
  numA = numA_;
  cout << "simSoft::setLJMinusPlusParams::numA: " << numA << endl;
}

void simSoft::setLJWCAparams(double LJcutoff_, long numA_) {
  setLJcutoff(LJcutoff_);
  numA = numA_;
  cout << "simSoft::setLJWCAparams::numA: " << numA << << endl;
}

void simSoft::addSelfPropulsion() {
  long f_nDim(nDim);
  double f_driving(driving);
  // use this to be able to set taup to 0
  double amplitude = 0.;
  if(taup != 0) {
    amplitude = sqrt(2.0 * dt / taup);
  }
  double *pAngle = angle.data();
  double *pForce = force.data();
	if(nDim == 2) {
    wrappedGaussNum generateWrappedNumbers(-PI, PI);
    std::for_each(randAngle.begin(), randAngle.end(), [&](double& val) {val = generateWrappedNumbers(0);});
    double *pRandAngle = randAngle.data();

    auto updateActiveNoise2D = [&](long pId) {
      pAngle[pId] += pRandAngle[pId];
      pAngle[pId] = pAngle[pId] + PI;
      pAngle[pId] = pAngle[pId] - 2.0 * PI * floor(pAngle[pId] / (2.0 * PI));
      pAngle[pId] = pAngle[pId] - PI;
      #pragma unroll (MAXDIM)
      for (long dim = 0; dim < f_nDim; dim++) {
        pForce[pId * f_nDim + dim] += f_driving * ((1 - dim) * cos(pAngle[pId]) + dim * sin(pAngle[pId]));
      }
    };

    std::for_each(angle.begin(), angle.end(), [&, i = 0](double&) mutable { updateActiveNoise2D(i++); });


  } else if(nDim == 3) {
    randNum generateGaussNumbers(-PI, PI);
    std::for_each(randAngle.begin(), randAngle.end(), [&](double& val) {val = generateGaussNumbers(0);});
    double *pRandAngle = randAngle.data();

    auto normalizeVector = [&](long pId) {
      auto norm = 0.0;
      #pragma unroll (MAXDIM)
      for (long dim = 0; dim < f_nDim; dim++) {
        norm += pRandAngle[pId * f_nDim + dim] * pRandAngle[pId * f_nDim + dim];
      }
      norm = sqrt(norm);
      #pragma unroll (MAXDIM)
      for (long dim = 0; dim < f_nDim; dim++) {
        pRandAngle[pId * f_nDim + dim] /= norm;
      }
    };

    std::for_each(randAngle.begin(), randAngle.end(), [&, i = 0](double&) mutable { normalizeVector(i++); });

    auto updateActiveNoise3D = [&](long pId) {
      pAngle[pId * f_nDim] += amplitude * (pAngle[pId * f_nDim + 1] * pRandAngle[pId * f_nDim + 2] - pAngle[pId * s_nDim + 2] * pRandAngle[pId * f_nDim + 1]);
      pAngle[pId * f_nDim + 1] += amplitude * (pAngle[pId * f_nDim + 2] * pRandAngle[pId * f_nDim] - pAngle[pId * s_nDim] * pRandAngle[pId * f_nDim + 2]);
      pAngle[pId * f_nDim + 2] += amplitude * (pAngle[pId * f_nDim] * pRandAngle[pId * f_nDim + 1] - pAngle[pId * s_nDim + 1] * pRandAngle[pId * f_nDim]);
      #pragma unroll (MAXDIM)
      for (long dim = 0; dim < f_nDim; dim++) {
        pForce[pId * f_nDim + dim] += f_driving * pAngle[pId * f_nDim + dim];
      }
    };

    std::for_each(angle.begin(), angle.end(), [&, i = 0](double&) mutable { updateActiveNoise3D(i++); });
  }
}

void simSoft::addVicsekAlignment() {
	if(nDim == 2) {
    int f_nDim(nDim);
    double f_dt(dt);
    double f_driving(driving);
    // use this to be able to set taup to 0
    double amplitude = 0.;
    if(taup != 0) {
      amplitude = sqrt(2.0 * dt / taup);
    }
    wrappedGaussNum generateWrappedNumbers(-PI, PI);
    std::for_each(randAngle.begin(), randAngle.end(), [&](double& val) {val = generateWrappedNumbers(0);});
    const double *pAlpha = alpha.data();
    double *randAngle = randAngle.data();
    double *pAngle = angle.data();
    double *pForce = force.data();

    auto updateVicsekAlignment2D = [&](long pId) {
      // overdamped equation for the angle with vicsek alignment as torque
      pAngle[pId] += pRandAngle[pId] + s_dt * pAlpha[pId];
      pAngle[pId] = pAngle[pId] + PI;
      pAngle[pId] = pAngle[pId] - 2.0 * PI * floor(pAngle[pId] / (2.0 * PI));
      pAngle[pId] = pAngle[pId] - PI;
      #pragma unroll (MAXDIM)
      for (long dim = 0; dim < f_nDim; dim++) {
        pForce[pId * f_nDim + dim] += f_driving * ((1 - dim) * cos(pAngle[pId]) + dim * sin(pAngle[pId]));
      }
    };

    std::for_each(r, r + numParticles, updateVicsekAlignment2D);
  } else {
    cout << "simSoft::addVicsekAlignment: only 2D vicsek alignment is implemented!" << endl;
  }
}

void simSoft::calcVicsekAlignment() {
  checkVicsekNeighbors();
  const double *pAngle = particleAngle.data();
  double *pAlpha = particleAlpha.data();
  if(nDim == 2) {
    auto calcVicsekAlignment2D = [&](long pId) {
    };
  } else {
    cout << "simSoft::calcVicsekAlignment: only 2D vicsek alignment is implemented!" << endl;
  }
}

void simSoft::calcInteraction() {
	const double *pRad = rad.data();
	const double *pPos = pos.data();
	double *pForce = force.data();
	double *pEnergy = energy.data();
  //cout << "dimGrid, dimBlock: " << dimGrid << ", " << dimBlock << endl;
  switch (simControl.neighborType) {
    case simControlStruct::neighborEnum::neighbor:
    kernelCalcParticleInteraction<<<dimGrid, dimBlock>>>(pRad, pPos, pForce, pEnergy);
    break;
    case simControlStruct::neighborEnum::allToAll:
    kernelCalcAllToAllParticleInteraction<<<dimGrid, dimBlock>>>(pRad, pPos, pForce, pEnergy);
    break;
    default:
    break;
  }
}

void simSoft::addFixedWallInteraction() {
  const double *pRad = rad.data();
	const double *pPos = pos.data();
	double *pForce = force.data();
	double *pEnergy = energy.data();
  kernelCalcFixedWallInteraction<<<dimGrid, dimBlock>>>(pRad, pPos, pForce, pEnergy);
}

void simSoft::calcForceEnergy() {
  calcParticleInteraction();
  switch (simControl.particleType) {
    case simControlStruct::particleEnum::active:
    addSelfPropulsion();
    break;
    case simControlStruct::particleEnum::vicsek:
    calcVicsekAlignment();
    addVicsekAlignment();
    break;
    default:
    break;
  }
  switch (simControl.boundaryType) {
    case simControlStruct::boundaryEnum::fixed:
    addFixedWallInteraction();
    break;
    default:
    break;
  }
}

double simSoft::getMaxUnbalancedForce() {
  std::vector<double> forceSquared(force.size());
  std::transform(force.begin(), force.end(), forceSquared.begin(), square());
  double maxUnbalancedForce = sqrt(std::reduce(forceSquared.begin(), forceSquared.end(), double(-1), std::maximum<double>()));
  forceSquared.clear();
  return maxUnbalancedForce;
}

double simSoft::getPotentialEnergy() {
  return std::reduce(energy.begin(), energy.end(), double(0), std::plus<double>());
}

double simSoft::getKineticEnergy() {
  std::transform(vel.begin(), vel.end(), squaredVel.begin(), square());
  return 0.5 * std::reduce(squaredVel.begin(), squaredVel.end(), double(0), std::plus<double>());
}

double simSoft::getTemperature() {
  return 2 * getKineticEnergy() / (nDim * numParticles);
}

double simSoft::getEnergy() {
  return (getPotentialEnergy() + getKineticEnergy());
}

void simSoft::adjustKineticEnergy(double prevEtot) {
  double scale, ekin = getKineticEnergy();
  double deltaEtot = getPotentialEnergy() + ekin;
  deltaEtot -= prevEtot;
  if(ekin > deltaEtot) {
    scale = sqrt((ekin - deltaEtot) / ekin);
    //cout << "deltaEtot: " << deltaEtot << " ekin - deltaEtot: " << ekin - deltaEtot << " scale: " << scale << endl;
    long f_nDim(nDim);
    double* pVel = vel.data();

    auto adjustVel = [&](long pId) {
      #pragma unroll (MAXDIM)
      for (long dim = 0; dim < f_nDim; dim++) {
        pVel[pId * f_nDim + dim] *= scale;
      }
    };

    cout << "simSoft::adjustKineticEnergy:: scale: " << scale << endl;
    std::for_each(r, r + numParticles, adjustVel);
  } else {
    cout << "simSoft::adjustKineticEnergy:: kinetic energy is less then change in total energy - no adjustment is made" << endl;
  }
}

//************************* neighbor update functions ***************************//
double simSoft::setDisplacementCutoff(double cutoff_) {
  switch (simControl.potentialType) {
    case simControlStruct::potentialEnum::harmonic:
    cutDistance = 1;
    break;
    case simControlStruct::potentialEnum::lennardJones:
    cutDistance = LJcutoff;
    break;
    case simControlStruct::potentialEnum::WCA:
    cutDistance = WCAcut;
    break;
    case simControlStruct::potentialEnum::doubleLJ:
    cutDistance = LJcutoff;
    break;
    case simControlStruct::potentialEnum::LJMinusPlus:
    cutDistance = LJcutoff;
    break;
    case simControlStruct::potentialEnum::LJWCA:
    cutDistance = LJcutoff;
    break;
    default:
    break;
  }
  cutDistance += cutoff_; // adimensional because it is used for the overlap (gap) between two particles
  cutoff = cutoff_ * getMeanParticleSize();
  cout << "simSoft::setDisplacementCutoff - cutDistance: " << cutDistance << " cutoff: " << cutoff << endl;
  return cutDistance;
}

void simSoft::resetUpdateCount() {
  updateCount = double(0);
  //cout << "simSoft::resetUpdateCount - updatCount " << updateCount << endl;
}

long simSoft::getUpdateCount() {
  return updateCount;
}

// this function is called after particleDisplacement has been computed
void simSoft::removeCOMDrift() {
  // compute drift on x
  std::vector<double> disp_x(disp.size() / 2);
  std::vector<long> idx(disp.size() / 2);
  std::iota(idx.begin(), idx.end(), 0);  // Generate sequence 0, 1, 2, ..., disp.size()/2 - 1
  for (size_t i = 0; i < idx.size(); ++i) {
      disp_x[i] = disp[idx[i] * 2];  // Manually gather every second element from disp
  }
  double drift_x = std::reduce(disp_x.begin(), disp_x.end(), double(0), std::plus<double>()) / numParticles;
  // compute drift on y
  std::vector<double> disp_y(disp.size() / 2);
  std::vector<long> idy(disp.size() / 2);
  for (size_t i = 0; i < idx.size(); ++i) {
    disp_x[i] = disp[idx[i] * 2 + 1];  // Manually gather every second element from disp
  }
  double drift_y = std::reduce(disp_y.begin(), disp_y.end(), double(0), std::plus<double>()) / numParticles;

  // subtract drift from current positions
  long f_nDim(nDim);
	double* pPos = pos.data();
	double* pBoxSize = boxSize.data();

  auto removeDrift = [&](long pId) {
		pPos[pId * f_nDim] -= drift_x;
    pPos[pId * f_nDim + 1] -= drift_y;
  };
  std::for_each(r, r + numParticles, removeDrift);
}

void simSoft::checkParticleDisplacement() {
  const double *pPos = pos.data();
  const double *pLastPos = lastPos.data();
  int *flag = dispFlag.data();
  kernelCheckParticleDisplacement<<<dimGrid,dimBlock>>>(pPos, pLastPos, flag, cutoff);
  int sumFlag = std::reduce(d_flag.begin(), d_flag.end(), int(0), std::plus<int>());
  if(sumFlag != 0) {
    calcParticleNeighborList(cutDistance);
    resetLastPositions();
    if(shift == true) {
      removeCOMDrift();
    }
    updateCount += 1;
  }
}

void simSoft::checkParticleNeighbors() {
  switch (simControl.neighborType) {
    case simControlStruct::neighborEnum::neighbor:
    checkParticleDisplacement();
    break;
    case simControlStruct::neighborEnum::allToAll:
    break;
    default:
    break;
  }
}

void simSoft::checkVicsekNeighbors() {
  const double *pPos = pos.data();
  const double *vLastPos = vicsekLastPos.data();
  int *vicsekFlag = vicsekFlag.data();
  kernelCheckParticleDisplacement<<<dimGrid,dimBlock>>>(pPos, vLastPos, vicsekFlag, Rvicsek);
  int sumFlag = std::reduce(d_vicsekFlag.begin(), d_vicsekFlag.end(), int(0), std::plus<int>());
  if(sumFlag != 0) {
    calcVicsekNeighborList();
    resetVicsekLastPositions();
  }
}

std::vector<long> simSoft::getNeighbors() {
  return neighborList;
}

std::vector<long> simSoft::getVicsekNeighbors() {
  return vicsekNeighborList;
}

void simSoft::calcNeighbors(double cutDistance) {
  switch (simControl.neighborType) {
    case simControlStruct::neighborEnum::neighbor:
    calcNeighborList(cutDistance);
    break;
    case simControlStruct::neighborEnum::allToAll:
    break;
    default:
    break;
  }
}

void simSoft::calcNeighborList(double cutDistance) {
  std::fill(d_partMaxNeighborList.begin(), d_partMaxNeighborList.end(), 0);
	std::fill(d_partNeighborList.begin(), d_partNeighborList.end(), -1L);
  syncParticleNeighborsToDevice();
  const double *pPos = pos.data();
	const double *pRad = rad.data();

  kernelCalcParticleNeighborList<<<dimGrid, dimBlock>>>(pPos, pRad, cutDistance);
  // compute maximum number of neighbors per particle
  if(cudaGetLastError()) cout << "simSoft::calcParticleNeighborList():: cudaGetLastError(): " << cudaGetLastError() << endl;
  partMaxNeighbors = std::reduce(d_partMaxNeighborList.begin(), d_partMaxNeighborList.end(), -1L, std::maximum<long>());
  syncParticleNeighborsToDevice();
  //cout << "simSoft::calcParticleNeighborList: maxNeighbors: " << partMaxNeighbors << endl;

  // if the neighbors don't fit, resize the neighbor list
  if ( partMaxNeighbors > partNeighborListSize ) {
		partNeighborListSize = pow(2, ceil(std::log2(partMaxNeighbors)));
    //cout << "simSoft::calcParticleNeighborList: neighborListSize: " << partNeighborListSize << endl;
		//Now create the actual storage and then put the neighbors in it.
		d_partNeighborList.resize(numParticles * partNeighborListSize);
		//Pre-fill the neighborList with -1
		std::fill(d_partNeighborList.begin(), d_partNeighborList.end(), -1L);
		syncParticleNeighborsToDevice();
		kernelCalcParticleNeighborList<<<dimGrid, dimBlock>>>(pPos, pRad, cutDistance);
	}
}

void simSoft::calcVicsekNeighborList() {
  std::fill(d_vicsekMaxNeighborList.begin(), d_vicsekMaxNeighborList.end(), 0);
	std::fill(d_vicsekNeighborList.begin(), d_vicsekNeighborList.end(), -1L);
  syncVicsekNeighborsToDevice();
  const double *pPos = pos.data();
	const double *pRad = rad.data();

  kernelCalcVicsekNeighborList<<<dimGrid, dimBlock>>>(pPos, pRad, Rvicsek);
  // compute maximum number of neighbors per particle
  if(cudaGetLastError()) cout << "simSoft::calcVicsekNeighborList():: cudaGetLastError(): " << cudaGetLastError() << endl;
  vicsekMaxNeighbors = std::reduce(d_vicsekMaxNeighborList.begin(), d_vicsekMaxNeighborList.end(), -1L, std::maximum<long>());
  syncVicsekNeighborsToDevice();
  //cout << "simSoft::calcVicsekNeighborList: vicsekMaxNeighbors: " << vicsekMaxNeighbors << endl;

  // if the neighbors don't fit, resize the neighbor list
  if ( vicsekMaxNeighbors > vicsekNeighborListSize ) {
		vicsekNeighborListSize = pow(2, ceil(std::log2(vicsekMaxNeighbors)));
    //cout << "simSoft::calcVicsekNeighborList: vicsekNeighborListSize: " << vicsekNeighborListSize << endl;
		//Now create the actual storage and then put the neighbors in it.
		d_vicsekNeighborList.resize(numParticles * vicsekNeighborListSize);
		//Pre-fill the neighborList with -1
		std::fill(d_vicsekNeighborList.begin(), d_vicsekNeighborList.end(), -1L);
		syncVicsekNeighborsToDevice();
		kernelCalcVicsekNeighborList<<<dimGrid, dimBlock>>>(pPos, pRad, Rvicsek);
	}
}

//************************** minimizer functions *****************************//
void simSoft::initFIRE(std::vector<double> &FIREparams, long minStep_, long numStep_, long numDOF_) {
  this->fire_ = new FIRE(this);
  if(FIREparams.size() == 7) {
    double a_start_ = FIREparams[0];
    double f_dec_ = FIREparams[1];
    double f_inc_ = FIREparams[2];
    double f_a_ = FIREparams[3];
    double fire_dt_ = FIREparams[4];
    double fire_dt_max_ = FIREparams[5];
    double a_ = FIREparams[6];
    this->fire_->initMinimizer(a_start_, f_dec_, f_inc_, f_a_, fire_dt_, fire_dt_max_, a_, minStep_, numStep_, numDOF_);
  } else {
    cout << "simSoft::initFIRE: wrong number of FIRE parameters, must be 7" << endl;
  }
  resetLastPositions();
}

void simSoft::setParticleMassFIRE() {
  this->fire_->setMass();
}

void simSoft::setTimeStepFIRE(double timeStep_) {
  this->fire_->setFIRETimeStep(timeStep_);
}

void simSoft::FIRELoop() {
  this->fire_->minimizerLoop();
}

//***************************** Langevin integrators ******************************//
void simSoft::initSoftLangevin(double Temp, double gamma, bool readState) {
  this->sim_ = new SoftLangevin(this, SimConfig(Temp));
  this->sim_->gamma = gamma;
  this->sim_->noise = sqrt(2. * Temp * gamma / dt);
  this->sim_->rand.resize(numParticles * nDim);
  std::fill(this->sim_->rand.begin(), this->sim_->rand.end(), double(0));
  resetLastPositions();
  setInitialPositions();
  if(readState == false) {
    this->sim_->injectKineticEnergy();
  }
  cout << "simSoft::initSoftLangevin:: current temperature: " << setprecision(12) << getTemperature() << endl;
}

void simSoft::softLangevinLoop(bool conserve) {
  this->sim_->integrate();
  if(conserve == true) {
    this->sim_->conserveMomentum();
  }
}

//***************************** NVE integrator *******************************//
void simSoft::initSoftNVE(double Temp, bool readState) {
  this->sim_ = new SoftNVE(this, SimConfig(Temp));
  resetLastPositions();
  setInitialPositions();
  shift = true;
  if(readState == false) {
    this->sim_->injectKineticEnergy();
  }
  cout << "simSoft::initSoftNVE:: current temperature: " << setprecision(12) << getTemperature() << endl;
}

void simSoft::softNVELoop() {
  this->sim_->integrate();
}

void simSoft::rescaleVelocities(double Temp) {
  double scale = sqrt(Temp / getTemperature());
  long f_nDim(nDim);
	double* pVel = vel.data();

  auto rescaleVel = [&](long pId) {
    #pragma unroll (MAXDIM)
		for (long dim = 0; dim < f_nDim; dim++) {
      pVel[pId * f_nDim + dim] *= scale;
    }
  };

  std::for_each(r, r + numParticles, rescaleVel);
}

//************************* Nose-Hoover integrator ***************************//
void simSoft::getNoseHooverParams(double &mass, double &damping) {
  mass = this->sim_->mass;
  damping = this->sim_->gamma;
  //cout << "SP2D::getNoseHooverParams:: damping: " << this->sim_->gamma << endl;
}

void simSoft::initSoftNoseHoover(double Temp, double mass, double gamma, bool readState) {
  this->sim_ = new SoftNoseHoover(this, SimConfig(Temp));
  this->sim_->mass = mass;
  this->sim_->gamma = gamma;
  resetLastPositions();
  setInitialPositions();
  shift = true;
  if(readState == false) {
    this->sim_->injectKineticEnergy();
  }
  cout << "simSoft::initSoftNoseHoover:: current temperature: " << setprecision(12) << getTemperature();
  cout << " mass: " << this->sim_->mass << ", damping: " << this->sim_->gamma << endl;
}

void simSoft::softNoseHooverLoop() {
  this->sim_->integrate();
}

//********************** functions for testing interaction *********************//
void simSoft::setTwoParticleTestPacking(double sigma0, double sigma1, double lx, double ly, double y0, double y1, double vel1) {
  // set particle radii
  rad[0] = 0.5 * sigma0;
  rad[1] = 0.5 * sigma1;
  boxSize[0] = lx;
  boxSize[1] = ly;
  // assign positions
  for (int pId = 0; pId < numParticles; pId++) {
    pos[pId * nDim] = lx * 0.5;
  }
  pos[0 * nDim + 1] = ly * y0;
  pos[1 * nDim + 1] = ly * y1;
  // assign velocity
  vel[1 * nDim + 1] = vel1;
}

void simSoft::setThreeParticleTestPacking(double sigma02, double sigma1, double lx, double ly, double y02, double y1, double vel1) {
  // set particle radii
  rad[0] = 0.5 * sigma02;
  rad[1] = 0.5 * sigma1;
  rad[2] = 0.5 * sigma02;
  boxSize[0] = lx;
  boxSize[1] = ly;
  // assign positions
  pos[0 * nDim] = lx * 0.35;
  pos[0 * nDim + 1] = ly * y02;
  pos[2 * nDim] = lx * 0.65;
  pos[2 * nDim + 1] = ly * y02;
  pos[1 * nDim] = lx * 0.5;
  pos[1 * nDim + 1] = ly * y1;
  // assign velocity
  vel[1 * nDim + 1] = vel1;
}

void simSoft::firstUpdate(double timeStep) {
  int f_nDim(nDim);
  double f_dt(timeStep);
	double* pPos = pos.data();
	double* pVel = vel.data();
	const double* pForce = force.data();

  auto firstUpdate = [&](long pId) {
    #pragma unroll (MAXDIM)
		for (long dim = 0; dim < f_nDim; dim++) {
      pVel[pId * f_nDim + dim] += 0.5 * f_dt * pForce[pId * f_nDim + dim];
      pPos[pId * f_nDim + dim] += f_dt * pVel[pId * f_nDim + dim];
    }
  };

  std::for_each(r, r + numParticles, firstUpdate);
}

void simSoft::secondUpdate(double timeStep) {
  int f_nDim(nDim);
  double f_dt(timeStep);
	double* pVel = vel.data();
	const double* pForce = force.data();

  auto firstUpdate = [&](long pId) {
    #pragma unroll (MAXDIM)
		for (long dim = 0; dim < f_nDim; dim++) {
      pVel[pId * f_nDim + dim] += 0.5 * f_dt * pForce[pId * f_nDim + dim];
    }
  };

  std::for_each(r, r + numParticles, firstUpdate);
}

void simSoft::testInteraction(double timeStep) {
  firstUpdate(timeStep);
  checkParticleNeighbors();
  calcParticleForceEnergy();
  secondUpdate(timeStep);
}