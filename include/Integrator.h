//
// Author: Francesco Arceri
// Date:   March 28 2025
//
// HEADER FILE FOR INTEGRATOR CLASS

#ifndef INTEGRATOR_H_
#define INTEGRATOR_H_

#include "simSoft.h"
#include <vector>

class simSoft;

class Integrator
{
public:
  simSoft * sp_;  //Pointer to the enclosing class

  Integrator() = default;
  Integrator(simSoft * spPtr);
  ~Integrator();
  // Global Integrator params
  double temp = 0.1;
  double mass = 1.;
  double gamma = 1.;
  double noise = 0.;
  std::vector<double> rand;

  //******************** functions for particle integrator ********************//
  void injectKineticEnergy();
  void updatePosition(double timeStep);
  void updateVelocity(double timeStep);
  void updateThermalVel();
  void conserveMomentum();
  // Integrator of Langevin dynamics
  void integrateLangevin();
  // Integrator at constant volume and energy
  void integrateNVE();
  // Functions for Nose-Hoover integrator
  void firstVelUpdateNH(double timeStep);
  void secondVelUpdateNH(double timeStep);
  void integrateNH();

};

#endif /* INTEGRATOR_H */
