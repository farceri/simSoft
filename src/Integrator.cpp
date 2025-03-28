//
// Author: Francesco Arceri
// Date:   March 28 2025
//
// DEFINITION OF INTEGRATION FUNCTIONS

#include "../include/Integrator.h"
#include "../include/defs.h"
#include "../include/simSoft.h"
#include <iostream>
#include <string>
#include <vector>
#include <cmath>
#include <random>

using namespace std;

//************************* soft particle langevin ***************************//
void SoftLangevin::integrate() {
  updateVelocity(0.5 * sp_->dt);
  updatePosition(sp_->dt);
  sp_->checkParticleNeighbors();
  sp_->calcParticleForceEnergy();
  updateThermalVel();
  updateVelocity(0.5 * sp_->dt);
}

void SoftLangevin::injectKineticEnergy() {
  double amplitude(sqrt(config.Tbath));
  // generate random numbers between 0 and noise for thermal noise
  thrust::counting_iterator<long> index_sequence_begin(lrand48());
  std::transform(index_sequence_begin, index_sequence_begin + sp_->numParticles * sp_->nDim, sp_->d_vel.begin(), gaussNum(0.f,amplitude));
}

void SoftLangevin::updateThermalVel() {
  // generate random numbers between 0 and 1 for thermal noise
  thrust::counting_iterator<long> index_sequence_begin(lrand48());
  std::transform(index_sequence_begin, index_sequence_begin + sp_->numParticles * sp_->nDim, rand.begin(), gaussNum(0.f,1.f));
  for (long pId = 0; pId < sp_->numParticles; pId++) {
		for (long dim = 0; dim < sp_->nDim; dim++) {
      sp_->force[pId * s_nDim + dim] += noise * rand[pId * sp_->nDim + dim] - gamma * pVel[pId * sp_->nDim + dim];
    }
  }
}

void SoftLangevin::updateVelocity(double timeStep) {
  for (long pId = 0; pId < sp_->numParticles; pId++) {
    for (long dim = 0; dim < sp_->nDim; dim++) {
      sp_->vel[pId * sp_->nDim + dim] += timeStep * sp_->force[pId * sp_->nDim + dim] / mass;
    }
  }
}

void SoftLangevin::updatePosition(double timeStep) {
  for (long pId = 0; pId < sp_->numParticles; pId++) {
    for (long dim = 0; dim < sp_->nDim; dim++) {
      sp_->pos[pId * sp_->nDim + dim] += timeStep * sp_->vel[pId * sp_->nDim + dim];
    }
  }
}

void SoftLangevin::conserveMomentum() {
  typedef std::vector<double>::iterator Iterator;
  strided_range<Iterator> vel_x(sp_->vel.begin(), sp_->vel.end(), 2);
  strided_range<Iterator> vel_y(sp_->vel.begin() + 1, sp_->vel.end(), 2);
  double meanVelx = std::reduce(vel_x.begin(), vel_x.end(), double(0), std::plus<double>()) / sp_->numParticles;
  double meanVely = std::reduce(vel_y.begin(), vel_y.end(), double(0), std::plus<double>()) / sp_->numParticles;

  for (long pId = 0; pId < sp_->numParticles; pId++) {
    for (long dim = 0; dim < sp_->nDim; dim++) {
      sp_->vel[pId * sp_->nDim + dim] -= ((1 - dim) * meanVelx + dim * meanVely);
    }
  }
}

//**************************** soft particle nve *****************************//
void SoftNVE::integrate() {
  updateVelocity(0.5 * sp_->dt);
  updatePosition(sp_->dt);
  sp_->checkParticleNeighbors();
  sp_->calcParticleForceEnergy();
  updateVelocity(0.5 * sp_->dt);
}

//************************ soft particle Nose Hoover **************************//
void SoftNoseHoover::integrate() {
  updateVelocity(0.5 * sp_->dt);
  updatePosition(sp_->dt);
  sp_->checkNeighbors();
  sp_->calcForceEnergy();
  updateThermalVel();
}

void SoftNoseHoover::updateVelocity(double timeStep) {
  // update nose hoover damping
  gamma += (timeStep / (2 * mass)) * (sp_->getKineticEnergy() - (sp_->nDim * sp_->numParticles + 1) * config.Tbath / 2);
  for (long pId = 0; pId < sp_->numParticles; pId++) {
    for (long dim = 0; dim < sp_->nDim; dim++) {
      sp_->vel[pId * sp_->nDim + dim] += timeStep * (sp_->force[pId * sp_->nDim + dim] - sp_->vel[pId * sp_->nDim + dim] * gamma);
    }
  }
}

void SoftNoseHoover::updateThermalVel() {
  double timeStep = 0.5 * sp_->dt;
  // update nose hoover damping
  gamma += (sp_->dt / (2 * mass)) * (sp_->getKineticEnergy() - (sp_->nDim * sp_->numParticles + 1) * config.Tbath / 2);
  for (long pId = 0; pId < sp_->numParticles; pId++) {
    for (long dim = 0; dim < sp_->nDim; dim++) {
      sp_->vel[pId * sp_->nDim + dim] = (sp_->vel[pId * sp_->nDim + dim] + timeStep * sp_->force[pId * sp_->nDim + dim]) / (1 + timeStep * gamma);
    }
  }
}