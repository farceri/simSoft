//
// Author: Francesco Arceri
// Date:   March 31 2025
//
// Testing simulation for pair interactions between particles.
// Tests are run for 2 and 3 particles with different potentials.
//
// Include C++ header files

#include "include/simSoft.h"
#include "include/FileIO.h"
#include "include/Integrator.h"
#include "include/defs.h"
#include <vector>
#include <string>
#include <iostream>
#include <iomanip>
#include <math.h>
#include <functional>
#include <utility>
#include <experimental/filesystem>

using namespace std;
using std::cout;

int main(int argc, char **argv) {
  // boolean variables for simulation setup
  bool read = false, readState = false, saveFinal = false, linSave = true;
  // boolean variables for potential type and force calculation
  bool doublelj = false, lj = false, wca = false, alltoall = true;
  // initialization variables
  long numParticles = atol(argv[2]), maxStep = atof(argv[3]), nDim = 2, num1 = 1;
  cout << "Number of particles: " << numParticles << endl;
  double timeStep = atof(argv[4]), ec = 1, ea = 1, eb = 1, eab = 1;
  double sigma0 = 1, sigma1 = 3, lx = 10, ly = 10, vel1 = -0.1, y0 = 0.2, y1 = 0.7;
  // other variables
  long step = 0, updateCount = 0, checkPointFreq = int(maxStep / 10);
  long linFreq = int(checkPointFreq / 10), saveEnergyFreq = int(linFreq / 10);
  double LJcut = 4, cutoff = 0.2, cutDistance, sigma, timeUnit, size;
  std::string inDir = argv[1], outDir, currentDir, dirSample, energyFile;
  // initialize sp object
  simSoft sp(numParticles, nDim);
  sp.setEnergyCostant(ec);

  // Potential type setting
  if(lj == true) {
    sp.setPotentialType(simControlStruct::potentialEnum::lennardJones);
    dirSample = "lj/";
    cout << "Setting Lennard-Jones potential" << endl;
    sp.setLJcutoff(LJcut);
  } else if(wca == true) {
    sp.setPotentialType(simControlStruct::potentialEnum::WCA);
    dirSample = "wca/";
    cout << "Setting WCA potential" << endl;
  } else if(doublelj == true) {
    sp.setPotentialType(simControlStruct::potentialEnum::doubleLJ);
    dirSample = "2lj/";
    cout << "Setting double Lennard-Jones potential" << endl;
    sp.setDoubleLJconstants(LJcut, ea, eab, eb, num1);
  } else {
    dirSample = "harmonic/";
    cout << "Setting Harmonic potential between particles and walls" << endl;
  }

  // Neighbor type setting
  if(alltoall == true) {
    sp.setNeighborType(simControlStruct::neighborEnum::allToAll);
  }

  // Input and output
  ioFile ioSP(&sp);
  if (read == true) {//keep running the same dynamics
    cout << "Read packing" << endl;
    inDir = inDir + dirSample;
    outDir = inDir;
    ioSP.readPackingFromDirectory(inDir, numParticles, nDim);
    if(readState == true) {
      ioSP.readState(inDir, numParticles, nDim);
    }
  } else {//start a new dyanmics
    cout << "Initialize new packing" << endl;
    if(numParticles == 3) {
      sp.setThreeParticleTestPacking(sigma0, sigma1, lx, ly, y0, y1, vel1);
      sp.printThreeParticles();
    } else if(numParticles == 2) {
      sp.setTwoParticleTestPacking(sigma0, sigma1, lx, ly, y0, y1, vel1);
      sp.printTwoParticles();
    } else {
      cout << "testInteraction only works for 2 or 3 particles!" << endl;
    }
    if(std::experimental::filesystem::exists(inDir + dirSample) == false) {
      std::experimental::filesystem::create_directory(inDir + dirSample);
    }
    outDir = inDir + dirSample;
  }
  std::experimental::filesystem::create_directory(outDir);
  // output file
  energyFile = outDir + "energy.dat";
  ioSP.openEnergyFile(energyFile);
  
  // initialization
  timeUnit = sigma0;//epsilon and mass are 1 sqrt(m sigma^2 / epsilon)
  timeStep = sp.setTimeStep(timeStep * timeUnit);
  cout << "Units - time: " << timeUnit << " space: " << sigma0 << endl;
  cout << "initial velocity on particle 1: " << vel1 << " time step: " << timeStep << endl;
  if(sp.getNeighborType() == simControlStruct::neighborEnum::neighbor) {
    sp.setDisplacementCutoff(cutoff);
    sp.calcNeighbors();
  }
  sp.calcForceEnergy();
  sp.resetUpdateCount();
  sp.setInitialPositions();

  // save initial configuration
  ioSP.savePacking(outDir);
  ioSP.saveNeighbors(outDir);

  // run integrator
  while(step != maxStep) {
    sp.verletLoop(timeStep);
    if(step % saveEnergyFreq == 0) {
      ioSP.saveEnergy(step, timeStep, numParticles);
      if(step % checkPointFreq == 0) {
        cout << "NVE: current step: " << step;
        cout << " K: " << sp.getKineticEnergy();
        cout << " U: " << sp.getPotentialEnergy();
        cout << " E: " << sp.getEnergy();
        if(sp.simControl.neighborType == simControlStruct::neighborEnum::neighbor) {
          updateCount = sp.getUpdateCount();
          if(step != 0 && updateCount > 0) {
            cout << " number of updates: " << updateCount << endl;
          } else {
            cout << " no updates" << endl;
          }
          sp.resetUpdateCount();
        } else {
          cout << endl;
        }
        if(saveFinal == true) {
          ioSP.savePacking(outDir);
          ioSP.saveNeighbors(outDir);
        }
      }
    }
    if(linSave == true) {
      if((step % linFreq) == 0) {
        currentDir = outDir + "/t" + std::to_string(step) + "/";
        std::experimental::filesystem::create_directory(currentDir);
        ioSP.saveState(currentDir);
        ioSP.saveForces(currentDir);
        ioSP.saveNeighbors(currentDir);
      }
    }
    step += 1;
  }
  
  // save final configuration
  if(saveFinal == true) {
    ioSP.savePacking(outDir);
    ioSP.saveNeighbors(outDir);
  }
  ioSP.closeEnergyFile();

  return 0;
}
