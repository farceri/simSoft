//
// Author: Francesco Arceri
// Date:   April 1 2025
//
// Running script for NVE dynamics
// Setup 1: run a sample for the first time readAndMakeNewDir = false, readAndSaveSameDir = false, runDynamics = false
// Setup 2: run and save in a new directory readAndMakeNewDir = true, readAndSaveSameDir = false, runDynamics = false, scaleVel = true (automatically set)
// Setup 3: run and save in the same directory readAndMakeNewDir = false, readAndSaveSameDir = true, runDynamics = false
// Setup 4: make a dynamics/ directory and run readAndMakeNewDir = false, readAndSaveSameDir = true, runDynamics = true
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

int main(int argc, char **argv) {
  // boolean variables for simulation setup
  bool readAndMakeNewDir = false, readAndSaveSameDir = false, runDynamics = false, scaleVel = false;
  bool alltoall = false, readState = true, saveFinal = true, logSave = false, linSave = true;
  // input variables
  long numParticles = atol(argv[2]), nDim = atol(argv[3]), initialStep = atof(argv[4]), maxStep = atof(argv[5]);
  double timeStep = atof(argv[6]), Tinject = atof(argv[7]); // Tinject is only used if scaleVel = true, use 0 otherwise
  // step variables
  long step = 0, updateCount = 0, firstDecade = 0, multiple = 1, saveFreq = 1;
  long checkPointFreq = int(maxStep / 10), linFreq = int(checkPointFreq / 10), saveEnergyFreq = int(linFreq / 10);
  // interaction variables
  double ec = 1, LJcut = 4, cutoff = 0.5, cutDistance, waveNum, sigma, timeUnit;
  std::string outDir, currentDir, dirSample, whichDynamics = "nve";
  std::string inDir = argv[1], potType = argv[8], energyFile;

  // initialize sp object
	simSoft sp(numParticles, nDim);
  sp.setEnergyCostant(ec);
  if(alltoall == true) {
    sp.setNeighborType(simControlStruct::neighborEnum::allToAll);
  }
  // set potential type
  if(potType == "lj") {
    whichDynamics = whichDynamics + argv[9] + "/";
    sp.setPotentialType(simControlStruct::potentialEnum::lennardJones);
    sp.setLJcutoff(LJcut);
  } else if(potType == "wca") {
    whichDynamics = whichDynamics + "/";
    sp.setPotentialType(simControlStruct::potentialEnum::WCA);
  } else {
    whichDynamics = whichDynamics + "/";
    cout << "Setting default harmonic potential" << endl;
  }
  // make new directory if needed
  if(std::experimental::filesystem::exists(inDir + whichDynamics) == false) {
    std::experimental::filesystem::create_directory(inDir + whichDynamics);
  }
  if(Tinject == 0) {
    dirSample = whichDynamics;
  } else {
    dirSample = whichDynamics + "T" + argv[7] + "/";
  }

  // set input and output
  ioFile ioSP(&sp);
  if (readAndSaveSameDir == true) {//keep running the same dynamics
    readState = true;
    inDir = inDir + dirSample;
    outDir = inDir;
    if(runDynamics == true) {
      if(logSave == true) {
        inDir = outDir;
        outDir = outDir + "dynamics-log/";
      }
      if(linSave == true) {
        inDir = outDir;
        outDir = outDir + "dynamics/";
      }
      if(std::experimental::filesystem::exists(outDir) == true) {
        inDir = outDir;
      } else {
        std::experimental::filesystem::create_directory(outDir);
      }
    }
  } else {//start a new dyanmics
    if(readAndMakeNewDir == true) {
      scaleVel = true; // automatically set velocity rescaling to keep new temperature constant
      readState = true;
      if(Tinject == 0) {
        outDir = inDir + "../" + dirSample;
      } else {
        outDir = inDir + "../../" + dirSample;
      }
    } else {
      if(std::experimental::filesystem::exists(inDir + whichDynamics) == false) {
        std::experimental::filesystem::create_directory(inDir + whichDynamics);
      }
      outDir = inDir + dirSample;
    }
    std::experimental::filesystem::create_directory(outDir);
  }
  cout << "inDir: " << inDir << endl << "outDir: " << outDir << endl;
  ioSP.readPackingFromDirectory(inDir, numParticles, nDim);
  if(readState == true) ioSP.readState(inDir, numParticles, nDim);
  // output file
  energyFile = outDir + "energy.dat";
  ioSP.openEnergyFile(energyFile);

  // initialization
  sigma = sp.getMeanParticleSize();
  timeUnit = sigma / sqrt(ec);//mass is 1 - sqrt(m sigma^2 / epsilon)
  timeStep = sp.setTimeStep(timeStep * timeUnit);
  cout << "Units - time: " << timeUnit << " space: " << sigma << " energy:" << ec << endl;
  cout << "Tinject: " << Tinject << " time step: " << timeStep << endl;
  if(scaleVel == true) {
    sp.initSoftNVERescale(Tinject);
  } else {
    sp.initSoftNVE(Tinject, readState);
  }
  sp.setDisplacementCutoff(cutoff);
  sp.calcNeighbors();
  sp.calcForceEnergy();
  sp.resetUpdateCount();
  waveNum = sp.getWaveNumber();

  // run integrator
  while(step != maxStep) {
    if(scaleVel == true) {
      sp.softNVERescaleLoop();
    } else {
      sp.softNVELoop();
    }
    if(step % saveEnergyFreq == 0) {
      ioSP.saveEnergyISF(step+initialStep, timeStep, numParticles, waveNum);
      if(step % checkPointFreq == 0) {
        cout << "NVE: current step: " << step;
        cout << " E/N: " << sp.getEnergy() / numParticles;
        cout << " T: " << sp.getTemperature();
        cout << " ISF: " << sp.getIsotropicISF(waveNum);
        if(sp.simControl.neighborType == simControlStruct::neighborEnum::neighbor) {
          updateCount = sp.getUpdateCount();
          if(step != 0 && updateCount > 0) {
            cout << " number of updates: " << updateCount << endl;
          } else {
            cout << " no updates" << endl;
          }
          sp.resetUpdateCount();
        } else cout << endl;
        if(saveFinal == true) {
          ioSP.savePacking(outDir);
          //ioSP.saveNeighbors(outDir);
        }
      }
    }
    if(logSave == true) {
      if(step > (multiple * checkPointFreq)) {
        saveFreq = 1;
        multiple += 1;
      }
      if((step - (multiple-1) * checkPointFreq) > saveFreq*10) {
        saveFreq *= 10;
      }
      if(((step - (multiple-1) * checkPointFreq) % saveFreq) == 0) {
        currentDir = outDir + "/t" + std::to_string(initialStep + step) + "/";
        std::experimental::filesystem::create_directory(currentDir);
        ioSP.saveState(currentDir);
        //ioSP.saveNeighbors(currentDir);
      }
    }
    if(linSave == true) {
      if((step % linFreq) == 0) {
        currentDir = outDir + "/t" + std::to_string(initialStep + step) + "/";
        std::experimental::filesystem::create_directory(currentDir);
        ioSP.saveState(currentDir);
        //ioSP.saveNeighbors(currentDir);
      }
    }
    step += 1;
  }

  // save final configuration
  if(saveFinal == true) {
    ioSP.savePacking(outDir);
    //ioSP.saveNeighbors(outDir);
  }
  ioSP.closeEnergyFile();

  return 0;
}
