//
// Author: Francesco Arceri
// Date:   March 28 2025
//
// FUNCTIONS FOR FIRE CLASS

#include "../include/FIRE.h"
#include "../include/simSoft.h"
#include <iostream>
#include <string>
#include <vector>
#include <cmath>
#include <random>

using namespace std;

//********************** constructor and deconstructor ***********************//
FIRE::FIRE(simSoft * spPtr):sp_(spPtr){
	// Note that mass is used only for particle-level FIRE
	mass.resize(sp_->numParticles * sp_->nDim);
	// Set variables to zero
	std::fill(mass.begin(), mass.end(), double(0));
}

FIRE::~FIRE() {
	mass.clear();
	velSquared.clear();
	forceSquared.clear();
};

// initilize the minimizer
void FIRE::initMinimizer(double a_start_, double f_dec_, double f_inc_, double f_a_, double fire_dt_, double fire_dt_max_, double a_, long minStep_, long numStep_, long numDOF_) {
	a_start = a_start_;
	f_dec = f_dec_;
	f_inc = f_inc_;
	f_a = f_a_;
	fire_dt = fire_dt_;
	fire_dt_max = fire_dt_max_;
	a = a_;
	minStep = minStep_;
	numStep = numStep_;
	velSquared.resize(numDOF_ * sp_->nDim);
	forceSquared.resize(numDOF_ * sp_->nDim);
}

//*************************** particle minimizer *****************************//
// update position and velocity in response to an applied force
void FIRE::updatePositionAndVelocity() {
	double totalForce(0);
	for (long pId = 0; pId < sp_->numParticles; pId++) {
		for (long dim = 0; dim < sp_->nDim; dim++) {
			sp_->vel[pId * sp_->nDim + dim] += 0.5 * fire_dt * sp_->force[pId * sp_->nDim + dim] / mass[pId * sp_->nDim + dim];
			sp_->pos[pId * sp_->nDim + dim] += fire_dt * sp_->vel[pId * sp_->nDim + dim];
			totalForce += sp_->force[pId * sp_->nDim + dim];
		}
	}
	if (totalForce == 0) {
		for (long pId = 0; pId < sp_->numParticles; pId++) {
			for (long dim = 0; dim < sp_->nDim; dim++) {
				sp_->vel[pId * sp_->nDim + dim] = 0;
			}
		}
	}
}

// update velocity in response to an applied force and return the maximum displacement in the previous step
void FIRE::updateVelocity() {
	double totalForce(0);
	for (long pId = 0; pId < sp_->numParticles; pId++) {
		for (long dim = 0; dim < sp_->nDim; dim++) {
			sp_->vel[pId * sp_->nDim + dim] += 0.5 * fire_dt * sp_->force[pId * sp_->nDim + dim] / mass[pId * sp_->nDim + dim];
			totalForce += sp_->force[pId * sp_->nDim + dim];
		}
		//If the total force on a particle is zero, then zero out the velocity as well
		if (totalForce == 0) {
			for (long pId = 0; pId < sp_->numParticles; pId++) {
				for (long dim = 0; dim < sp_->nDim; dim++) {
					sp_->vel[pId * sp_->nDim + dim] = 0;
				}
			}
		}
	}
}

// bend the velocity towards the force
void FIRE::bendVelocityTowardsForce() {
	double velNormSquared = 0., forceNormSquared = 0.;
	// get the dot product between the velocity and the force
	double vDotF = double(std::inner_product(sp_->d_vel.begin(), sp_->d_vel.end(), sp_->d_force.begin(), double(0)));
	//cout << "FIRE::bendVelocityTowardsForceFIRE: vDotF = " << setprecision(precision) << vDotF << endl;
	if (vDotF < 0) {
		// if vDotF is negative, then we are going uphill, so let's stop and reset
		std::fill(sp_->d_vel.begin(), sp_->d_vel.end(), double(0));
		numStep = 0;
		fire_dt = std::max(fire_dt * f_dec, fire_dt_max / 10); // go to a shorter dt
		a = a_start; // start fresh with a more radical mixing between force and velocity
	} else if (numStep > minStep) {
		// if enough time has passed then let's start to increase the inertia
		fire_dt = std::min(fire_dt * f_inc, fire_dt_max);
		a *= f_a; // increase the inertia
	}
	// calculate the ratio of the norm squared of the velocity and the force
  std::transform(sp_->d_vel.begin(), sp_->d_vel.end(), d_velSquared.begin(), square());
  std::transform(sp_->d_force.begin(), sp_->d_force.end(), d_forceSquared.begin(), square());
	velNormSquared = std::reduce(velSquared.begin(), velSquared.end(), double(0), std::plus<double>());
	forceNormSquared = std::reduce(forceSquared.begin(), forceSquared.end(), double(0), std::plus<double>());
	// check FIRE convergence
	if (forceNormSquared == 0) {
		// if the forceNormSq is zero, then there is no force and we are done, so zero out the velocity
		cout << "FIRE::bendVelocityTowardsForceFIRE: forceNormSquared is zero" << endl;
		std::fill(sp_->d_vel.begin(), sp_->d_vel.end(), double(0));
	} else {
		double velForceNormRatio = sqrt(velNormSquared / forceNormSquared);
		double func_a(a);
		double* pVel = sp_->vel.data();
		const double* pForce = sp_->force.data();
		
		auto perDOFBendParticleVelocity = [&](long i) {
			pVel[i] = (1 - func_a) * pVel[i] + func_a * pForce[i] * velForceNormRatio;
		};

		std::for_each(sp_->vel.begin(), sp_->vel.end(), [&, i = 0](double&) mutable { perDOFBendParticleVelocity(i++); });
	}
}

// set the mass for each degree of freedom
void FIRE::setMass() {
	d_mass.resize(sp_->numParticles * sp_->nDim);
	for (long pId = 0; pId < sp_->numParticles; pId++) {
		for (long dim = 0; dim < sp_->nDim; dim++) {
			mass[pId * sp_->nDim + dim] = PI / (sp_->d_rad[pId] * sp_->d_rad[pId]);
		}
	}
}

// set FIRE time step
void FIRE::setFIRETimeStep(double fire_dt_) {
	fire_dt = fire_dt_;
	fire_dt_max = 10 * fire_dt_;
}

// Run the inner loop of the FIRE algorithm
void FIRE::minimizerLoop() {
	// Move the system forward, based on the previous velocities and forces
	updatePositionAndVelocity();
	// Calculate the new set of forces at the new step
  	sp_->checkNeighbors();
	sp_->calcForceEnergy();
	// update the velocity based on the current forces
	updateVelocity();
	// Bend the velocity towards the force
	bendVelocityTowardsForce();
	// Increase the number of steps since the last restart
	numStep++;
}
