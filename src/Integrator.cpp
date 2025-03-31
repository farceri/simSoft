//
// Author: Francesco Arceri
// Date:   March 28 2025
//
// DEFINITION OF INTEGRATION FUNCTIONS

#include "../include/Integrator.h"
#include "../include/defs.h"
#include <iostream>
#include <string>
#include <vector>
#include <cmath>
#include <random>

using namespace std;

//********************** constructor and deconstructor ***********************//
Integrator::Integrator(simSoft * spPtr):sp_(spPtr){
	// Vector for extracting white noise
	rand.resize(sp_->numParticles * sp_->nDim);
	// Set variables to zero
	std::fill(rand.begin(), rand.end(), double(0));
}

Integrator::~Integrator() {
	rand.clear();
};

void Integrator::injectKineticEnergy() {
  double amplitude(sqrt(temp));
  // generate random numbers between 0 and noise for thermal noise
  gaussNum generateGaussNumbers(0.f, 1.f);
  std::for_each(sp_->vel.begin(), sp_->vel.end(), [&](double& val) {val = generateGaussNumbers(0);});
}

void Integrator::updateThermalVel() {
  // generate random numbers between 0 and 1 for thermal noise
  gaussNum generateGaussNumbers(0.f, 1.f);
  std::for_each(rand.begin(), rand.end(), [&](double& val) {val = generateGaussNumbers(0);});
  for (long pId = 0; pId < sp_->numParticles; pId++) {
		for (long dim = 0; dim < sp_->nDim; dim++) {
      sp_->force[pId * sp_->nDim + dim] += noise * rand[pId * sp_->nDim + dim] - gamma * sp_->vel[pId * sp_->nDim + dim];
    }
  }
}

void Integrator::updateVelocity(double timeStep) {
  for (long pId = 0; pId < sp_->numParticles; pId++) {
    for (long dim = 0; dim < sp_->nDim; dim++) {
      sp_->vel[pId * sp_->nDim + dim] += timeStep * sp_->force[pId * sp_->nDim + dim] / mass;
    }
  }
}

void Integrator::updatePosition(double timeStep) {
  for (long pId = 0; pId < sp_->numParticles; pId++) {
    for (long dim = 0; dim < sp_->nDim; dim++) {
      sp_->pos[pId * sp_->nDim + dim] += timeStep * sp_->vel[pId * sp_->nDim + dim];
    }
  }
}

void Integrator::conserveMomentum() {
  std::vector<double> vel_x(sp_->vel.size() / 2);
  std::vector<double> vel_y(sp_->vel.size() / 2);
  std::vector<long> idx(sp_->vel.size() / 2);
  std::iota(idx.begin(), idx.end(), 0);  // Generate sequence 0, 1, 2, ..., vel.size()/2 - 1
  for (size_t i = 0; i < idx.size(); ++i) {
      vel_x[i] = sp_->vel[idx[i] * 2];  // Manually gather every even element from vel
      vel_y[i] = sp_->vel[idx[i] * 2 + 1];  // Manually gather every odd element from vel
  }
  double meanVelx = std::reduce(vel_x.begin(), vel_x.end(), double(0), std::plus<double>()) / sp_->numParticles;
  double meanVely = std::reduce(vel_y.begin(), vel_y.end(), double(0), std::plus<double>()) / sp_->numParticles;

  for (long pId = 0; pId < sp_->numParticles; pId++) {
    for (long dim = 0; dim < sp_->nDim; dim++) {
      sp_->vel[pId * sp_->nDim + dim] -= ((1 - dim) * meanVelx + dim * meanVely);
    }
  }
}

//************************* soft particle langevin ***************************//
void Integrator::integrateLangevin() {
  updateVelocity(0.5 * sp_->dt);
  updatePosition(sp_->dt);
  sp_->checkNeighbors();
  sp_->calcForceEnergy();
  updateThermalVel();
  updateVelocity(0.5 * sp_->dt);
}

//**************************** soft particle nve *****************************//
void Integrator::integrateNVE() {
  updateVelocity(0.5 * sp_->dt);
  updatePosition(sp_->dt);
  sp_->checkNeighbors();
  sp_->calcForceEnergy();
  updateVelocity(0.5 * sp_->dt);
}

void Integrator::firstVelUpdateNH(double timeStep) {
  // update nose hoover damping
  gamma += (timeStep / (2 * mass)) * (sp_->getKineticEnergy() - (sp_->nDim * sp_->numParticles + 1) * temp / 2);
  for (long pId = 0; pId < sp_->numParticles; pId++) {
    for (long dim = 0; dim < sp_->nDim; dim++) {
      sp_->vel[pId * sp_->nDim + dim] += timeStep * (sp_->force[pId * sp_->nDim + dim] - sp_->vel[pId * sp_->nDim + dim] * gamma);
    }
  }
}

void Integrator::secondVelUpdateNH(double timeStep) {
  // update nose hoover damping
  gamma += (sp_->dt / (2 * mass)) * (sp_->getKineticEnergy() - (sp_->nDim * sp_->numParticles + 1) * temp / 2);
  for (long pId = 0; pId < sp_->numParticles; pId++) {
    for (long dim = 0; dim < sp_->nDim; dim++) {
      sp_->vel[pId * sp_->nDim + dim] = (sp_->vel[pId * sp_->nDim + dim] + timeStep * sp_->force[pId * sp_->nDim + dim]) / (1 + timeStep * gamma);
    }
  }
}

//************************ soft particle Nose Hoover **************************//
void Integrator::integrateNH() {
  firstVelUpdateNH(0.5 * sp_->dt);
  updatePosition(sp_->dt);
  sp_->checkNeighbors();
  sp_->calcForceEnergy();
  secondVelUpdateNH(0.5 * sp_->dt);
}