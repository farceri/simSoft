//
// Author: Francesco Arceri
// Date:   March 30 2025
//
// Compression for producing particle samples (packings) in 2 and 3 dimensions.
// Samples are initialized with the FIRE minimizer in harmonic potential to suppress large forces.
// Each compression step consists of running dynamics for maxStep and then shrink the simulation box.
// Input directory is argv[1], es. mkdir samples/
//
// Sample run
// ./compress /path/to/directory/twoSpecies/samples/ 1024 2 1e-04 1e-03 4 1
//
// Include C++ header files

#include "include/simSoft.h"
#include "include/FileIO.h"
#include "include/Integrator.h"
#include "include/FIRE.h"
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

int main(int argc, char **argv) {
  // boolean variables for simulation setup
  bool read = false, readState = false, nve = true, noseHoover = false;
  // boolean variables for potential type and force calculation
  bool lj = false, wca = false, alltoall = false;
  // input variables
  long numParticles = atol(argv[2]), nDim = atol(argv[3]);
  double dt = atof(argv[4]), Tinject = atof(argv[5]), lx = atof(argv[6]), ly = atof(argv[7]), lz;
  if(nDim == 3) lz = atof(argv[8]);
  else lz = 0;
  // minimizer variables
  long iteration = 0, maxIterations = 1e05, minStep = 20, numStep = 0;
  double FIREStep = 1e-02, forceTollerance = 1e-08;
  // compression iteration variables
  long maxSearchStep = 1500, searchStep = 0, maxStep = 1e05, step = 0;
  long printFreq = int(maxStep / 10), updateCount = 0, saveEnergyFreq = int(printFreq / 10);
  double previousPhi, currentPhi, scaleFactor, prevEnergy = 0, deltaPhi = 1e-02, phi0 = 0.01, phiTh = 0.72;
  // other variables
  double LJcut = 4, poly = 0.2, ec = 1, cutoff = 0.5, mass = 1, damping = 1;
  if(nDim == 3) LJcut = 2.5;
  double size, cutDistance, timeStep, timeUnit, sigma;
  long num1 = int(numParticles / 2);
  std::string currentDir, outDir = argv[1], inDir, energyFile;
  std::vector<double> boxSize(nDim);
  // FIRE minimizer paramaters: a_start, f_dec, f_inc, f_a, dt, dt_max, a
  std::vector<double> particleFIREparams = {0.2, 0.5, 1.1, 0.99, FIREStep, 10*FIREStep, 0.2};

	// Initialize sp object
	simSoft sp(numParticles, nDim);
  sp.setEnergyCostant(ec);
  if(alltoall == true) {
    sp.setNeighborType(simControlStruct::neighborEnum::allToAll);
  }
  ioFile ioSP(&sp);
  std::experimental::filesystem::create_directory(outDir);
  // Read initial configuration: argv[9] is the name of the directory to read.
  // Directories are saved with the value of density (packing fraction) of the sytstem, ex. 0.60
  if(read == true) {
    inDir = argv[9];
    inDir = outDir + inDir + "/";
    ioSP.readPackingFromDirectory(inDir, numParticles, nDim);
    sigma = sp.getMeanParticleSize();
    if(readState == true) {
      ioSP.readState(inDir, numParticles, nDim);
    }
  } else { // Use harmonic potential for FIRE minimization
    // Initialize polydisperse packing
    sp.setPolyRandomParticles(phi0, poly, lx, ly, lz);
    //sp.setScaledMonoRandomParticles(phi0, lx, ly, lz);
    //sp.setScaledBiRandomParticles(phi0, lx, ly, lz);
    sp.scalePacking();
    sigma = sp.getMeanParticleSize();
    sp.initFIRE(particleFIREparams, minStep, numStep, numParticles);
    sp.setMassFIRE();
    sp.setDisplacementCutoff(cutoff);
    sp.calcNeighbors();
    sp.calcForceEnergy();
    sp.resetUpdateCount();
    sp.setInitialPositions();
    cout << "Generate initial configurations with FIRE \n" << endl;
    while((sp.getMaxUnbalancedForce() > forceTollerance) && (iteration != maxIterations)) {
      sp.FIRELoop();
      if(iteration % printFreq == 0 && iteration != 0) {
        cout << "FIRE: iteration: " << iteration;
        cout << " maxUnbalancedForce: " << setprecision(precision) << sp.getMaxUnbalancedForce();
        cout << " energy: " << sp.getEnergy() << "\n" << endl;
      }
      iteration += 1;
    }
    cout << "FIRE: last iteration: " << iteration;
    cout << " maxUnbalancedForce: " << setprecision(precision) << sp.getMaxUnbalancedForce();
    cout << " energy: " << setprecision(precision) << sp.getEnergy() << "\n" << endl;
  }

  // Set potential interaction for compression protocol
  if(lj == true) {
    sp.setPotentialType(simControlStruct::potentialEnum::lennardJones);
    cout << "Setting Lennard-Jones potential" << endl;
    sp.setLJcutoff(LJcut);
  } else if(wca == true) {
    sp.setPotentialType(simControlStruct::potentialEnum::WCA);
    cout << "Setting WCA potential" << endl;
  }

  // Save initial configuration
  currentDir = outDir + "initial/";
  std::experimental::filesystem::create_directory(currentDir);
  ioSP.savePacking(currentDir);

  // Quasistatic compression at finite temperature
  currentPhi = sp.getPackingFraction();
  cout << "Current density: " << currentPhi << ", average size: " << sigma << endl;
  previousPhi = currentPhi;
  timeUnit = sigma / sqrt(ec);
  timeStep = sp.setTimeStep(dt*timeUnit);
  if(nve == true) {
    cout << "NVE: Time step: " << timeStep << ", Tinject: " << Tinject << endl;
    sp.initSoftNVE(Tinject, readState);
  } else if(noseHoover == true) {
    cout << "Nose-Hoover: Time step: " << timeStep << ", Tinject: " << Tinject << endl;
    sp.initSoftNoseHoover(Tinject, mass, damping, readState);
  } else {
    cout << "Langevin: Time step: " << timeStep << ", Tinject: " << Tinject << ", damping: " << damping << endl;
    sp.initSoftLangevin(Tinject, damping, readState);
  }

  // Compression loop
  while (searchStep < maxSearchStep) {
    currentDir = outDir + std::to_string(sp.getPackingFraction()).substr(0,5) + "/";
    std::experimental::filesystem::create_directory(currentDir);
    energyFile = currentDir + "energy.dat";
    ioSP.openEnergyFile(energyFile);
    sp.setDisplacementCutoff(cutoff); // Need to do this again if the potential changed
    sp.calcNeighbors();
    sp.calcForceEnergy();
    sp.resetUpdateCount();
    sp.setInitialPositions();

    // Remove energy injected by compression for NVE
    if(searchStep != 0 && nve == true) {
      sp.calcNeighbors();
      sp.calcForceEnergy();
      cout << "Energy after compression - E/N: " << sp.getEnergy() / numParticles << endl;
      sp.adjustKineticEnergy(prevEnergy);
      sp.calcForceEnergy();
      cout << "Energy after adjustment - E/N: " << sp.getEnergy() / numParticles << endl;
    }

    // Run dynamics to equilibrate compressed configuration
    step = 0;
    while(step != maxStep) {
      if(nve == true) {
        sp.softNVELoop();
      } else if(noseHoover == true) {
        sp.softNoseHooverLoop();
      } else {
        sp.softLangevinLoop();
      }
      if(step % saveEnergyFreq == 0) {
        ioSP.saveEnergy(step, timeStep, numParticles);
      }
      if(step % printFreq == 0) {
        cout << "Compression: current step: " << step;
        cout << " E/N: " << sp.getEnergy() / numParticles;
        cout << " T: " << sp.getTemperature();
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
      }
      step += 1;
    }
    prevEnergy = sp.getEnergy();
    cout << "Energy before compression - E/N: " << prevEnergy / numParticles << endl;

    // Save equilibrated configuration
    ioSP.savePacking(currentDir);
    ioSP.saveNeighbors(currentDir);
    ioSP.closeEnergyFile();

    // Check if target density is met
    if(currentPhi > phiTh) {
      cout << "\nTarget density met, current density: " << currentPhi << endl;
      searchStep = maxSearchStep; // exit condition
    } else { // Compute scaleFactor to increase particle density
      if(nDim == 2) {
        scaleFactor = sqrt((currentPhi + deltaPhi) / currentPhi);
      } else if(nDim == 3) {
        scaleFactor = cbrt((currentPhi + deltaPhi) / currentPhi);
      } else {
        cout << "ScaleFactor: only dimensions 2 and 3 are allowed!" << endl;
      }
      sp.scaleParticles(scaleFactor);
      sp.scalePacking();
      currentPhi = sp.getPackingFraction();
      boxSize = sp.getBoxSize();
      if(nDim == 2) {
        cout << "\nNew phi: " << currentPhi << " Lx: " << boxSize[0] << " Ly: " << boxSize[1] << " scale: " << scaleFactor << endl;
      } else if(nDim == 3) {
        cout << "\nNew phi: " << currentPhi << " Lx: " << boxSize[0] << " Ly: " << boxSize[1] << " Lz: " << boxSize[2] << " scale: " << scaleFactor << endl;
      } else {
        cout << "Box rescaling: only dimensions 2 and 3 are allowed!" << endl;
      }
      searchStep += 1;
    }
  }
  return 0;
}
