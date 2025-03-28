//
// Author: Francesco Arceri
// Date:   March 28 2025
//
// HEADER FILE FOR INTEGRATOR CLASS

#ifndef INTEGRATOR_H
#define INTEGRATOR_H

#include "simSoft.h"
#include "defs.h"
#include <vector>

class Integrator;

class Integrator // initializer
{
public:
  double Tbath = 0.;
  SimConfig() = default;
  SimConfig(double Tbath):Tbath(Tbath){}
};

class SimInterface // integration functions
{
public:
  simSoft * sp_;
  SimConfig config;
  double mass = 1.;
  double gamma = 1.;
  double noise = 0.;

  SimInterface() = default;
  SimInterface(simSoft * spPtr, SimConfig config):sp_(spPtr),config(config){}
  virtual ~SimInterface() = default;

  virtual void injectKineticEnergy() = 0;
  virtual void updatePosition(double timeStep) = 0;
  virtual void updateVelocity(double timeStep) = 0;
  virtual void updateThermalVel() = 0;
  virtual void conserveMomentum() = 0;
  virtual void integrate() = 0;
};

//********************* integrators for soft particles ***********************//
// Soft particle Langevin integrator child of SimInterface
class SoftLangevin: public SimInterface
{
public:
  SoftLangevin() = default;
  SoftLangevin(SP2D * spPtr, SimConfig config) : SimInterface:: SimInterface(spPtr, config){;}

  virtual void injectKineticEnergy();
  virtual void updatePosition(double timeStep);
  virtual void updateVelocity(double timeStep);
  virtual void updateThermalVel();
  virtual void conserveMomentum();
  virtual void integrate();
};

// Soft particle NVE integrator child of SoftLangevin
class SoftNVE: public SoftLangevin
{
public:
  SoftNVE() = default;
  SoftNVE(SP2D * spPtr, SimConfig config) : SoftLangevin:: SoftLangevin(spPtr, config){;}

  virtual void integrate();
};

// Soft particle Nose-Hoover integrator child of SoftLangevin
class SoftNoseHoover: public SoftLangevin
{
public:
  SoftNoseHoover() = default;
  SoftNoseHoover(SP2D * spPtr, SimConfig config) : SoftLangevin:: SoftLangevin(spPtr, config){;}

  virtual void updateVelocity(double timeStep);
  virtual void updateThermalVel();
  virtual void integrate();
};

#endif // SIMULATOR_H //
