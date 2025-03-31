//
// Author: Francesco Arceri
// Date:   March 28 2025
//
// HEADER FILE FOR FIRE CLASS

#ifndef FIRE_H_
#define FIRE_H_

#include "simSoft.h"
#include <vector>

class simSoft;

class FIRE
{
public:
	simSoft * sp_;  //Pointer to the enclosing class

	FIRE() = default;
	FIRE(simSoft * spPtr);
	~FIRE();
	// Global FIRE params
	double a_start;
	double a;
	double f_dec;
	double f_inc;
	double f_a;
	double fire_dt_max;
	double fire_dt;
	double cutDistance;
	long minStep;
	long numStep;
	// FIRE variables
	std::vector<double> mass;
	std::vector<double> velSquared;
	std::vector<double> forceSquared;

	// initialize minimizer for particles
	void initMinimizer(double a_start_, double f_dec_, double f_inc_, double f_a_, double fire_dt_, 
		double dt_max_, double a_, long minStep_, long numStep_, long numDOF_);

	//******************** functions for particle minimizer ********************//
	// update position and velocity in response to an applied force
	void updatePositionAndVelocity();
	// update velocity in response to an applied force and return the maximum displacement in the previous step
	void updateVelocity();
	// bend the velocity towards the force
	void bendVelocityTowardsForce();
	// set the mass for each degree of freedom
	void setMass();
	// set FIRE time step
	void setFIRETimeStep(double timeStep);
	// fire minimizer loop for particles
	void minimizerLoop();

};

#endif /* FIRE_H_ */
