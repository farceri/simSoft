//
// Author: Francesco Arceri
// Date: March 28 2025
//

#ifndef FILEIO_H
#define FILEIO_H

#include "simSoft.h"
#include "defs.h"
#include <iostream>
#include <string>
#include <vector>
#include <fstream>
#include <cmath>

using namespace std;

class ioFile
{
public:
  ifstream inputFile;
  ofstream outputFile;
  ofstream energyFile;
  ofstream memoryFile;
  ofstream corrFile;
  simSoft * sp_;

  ioFile() = default;
  ioFile(simSoft * spPtr) {
    this->sp_ = spPtr;
  }

  // open file and check if it throws an error
  void openInputFile(string fileName) {
    inputFile = ifstream(fileName.c_str());
    if (!inputFile.is_open()) {
      cerr << "ioFile::openInputFile: error: could not open input file " << fileName << endl;
      exit(1);
    }
  }

  void openOutputFile(string fileName) {
    outputFile = ofstream(fileName.c_str());
    if (!outputFile.is_open()) {
      cerr << "ioFile::openOutputFile: error: could not open input file " << fileName << endl;
      exit(1);
    }
  }

  void openEnergyFile(string fileName) {
    energyFile = ofstream(fileName.c_str());
    if (!energyFile.is_open()) {
      cerr << "ioFile::openEnergyFile: error: could not open input file " << fileName << endl;
      exit(1);
    }
  }
  
  void closeEnergyFile() {
    energyFile.close();
  }

  void saveEnergy(long step, double timeStep, long numParticles) {
    double epot = sp_->getPotentialEnergy();
    double ekin = sp_->getKineticEnergy();
    double etot = epot + ekin;
    energyFile << step + 1 << "\t" << (step + 1) * timeStep << "\t";
    energyFile << setprecision(precision) << epot / numParticles << "\t";
    energyFile << setprecision(precision) << ekin / numParticles << "\t";
    energyFile << setprecision(precision) << etot / numParticles << endl;
  }

  //////////////////////////// write this function in array form ////////////////////////////
  double read0DFile(string fileName) {
    double data;
    this->openInputFile(fileName);
    string inputString;
    getline(inputFile, inputString);
    sscanf(inputString.c_str(), "%lf", &data);
    inputFile.close();
    return data;
  }

  void save0DFile(string fileName, double data) {
    this->openOutputFile(fileName);
    outputFile << setprecision(precision) << data << endl;
    outputFile.close();
  }

  std::vector<long> read1DIndexFile(string fileName, long numRows) {
    std::vector<long> data;
    this->openInputFile(fileName);
    string inputString;
    long tmp;
    for (long row = 0; row < numRows; row++) {
      getline(inputFile, inputString);
      sscanf(inputString.c_str(), "%ld", &tmp);
      data.push_back(tmp);
    }
    inputFile.close();
    return data;
  }

  void save1DIndexFile(string fileName, std::vector<long> data) {
    this->openOutputFile(fileName);
    long numRows = data.size();
    for (long row = 0; row < numRows; row++) {
      //sprintf(outputFile, "%ld \n", data[row]);
      outputFile << setprecision(precision) << data[row] << endl;
    }
    outputFile.close();
  }

  std::vector<double> read1DFile(string fileName, long numRows) {
    std::vector<double> data;
    this->openInputFile(fileName);
    string inputString;
    double tmp;
    for (long row = 0; row < numRows; row++) {
      getline(inputFile, inputString);
      sscanf(inputString.c_str(), "%lf", &tmp);
      data.push_back(tmp);
    }
    inputFile.close();
    return data;
  }

  void save1DFile(string fileName, std::vector<double> data) {
    this->openOutputFile(fileName);
    long numRows = data.size();
    for (long row = 0; row < numRows; row++) {
      //sprintf(outputFile, "%lf \n", data[row]);
      outputFile << setprecision(precision) << data[row] << endl;
    }
    outputFile.close();
  }

  std::vector<double> read2DFile(string fileName, long numRows) {
    std::vector<double> data;
    this->openInputFile(fileName);
    string inputString;
    double data1, data2;
    for (long row = 0; row < numRows; row++) {
      getline(inputFile, inputString);
      sscanf(inputString.c_str(), "%lf %lf", &data1, &data2);
      data.push_back(data1);
      data.push_back(data2);
    }
    inputFile.close();
    return data;
  }

  void save2DFile(string fileName, std::vector<double> data, long numCols) {
    this->openOutputFile(fileName);
    long numRows = int(data.size()/numCols);
    for (long row = 0; row < numRows; row++) {
      for(long col = 0; col < numCols; col++) {
        outputFile << setprecision(precision) << data[row * numCols + col] << "\t";
      }
      outputFile << endl;
    }
    outputFile.close();
  }

  void save2DIndexFile(string fileName, std::vector<long> data, long numCols) {
    this->openOutputFile(fileName);
    if(numCols == 0) {
      cout << "save2dIndexFile: empty array!" << endl;
      return;
    } else {
      long numRows = int(data.size()/numCols);
      for (long row = 0; row < numRows; row++) {
        for(long col = 0; col < numCols; col++) {
          outputFile << data[row * numCols + col] << "\t";
        }
        outputFile << endl;
      }
      outputFile.close();
    }
  }

  std::vector<double> read3DFile(string fileName, long numRows) {
    std::vector<double> data;
    this->openInputFile(fileName);
    string inputString;
    double data1, data2, data3;
    for (long row = 0; row < numRows; row++) {
      getline(inputFile, inputString);
      sscanf(inputString.c_str(), "%lf %lf %lf", &data1, &data2, &data3);
      data.push_back(data1);
      data.push_back(data2);
      data.push_back(data3);
    }
    inputFile.close();
    return data;
  }

  std::vector<double> readBoxSize(string dirName, long nDim_) {
    std::vector<double> boxSize_(nDim_);
    boxSize_ = read1DFile(dirName + "boxSize.dat", nDim_);
    cout << "FileIO::readBoxSize: " << " box-Lx: " << boxSize_[0] << ", Ly: " << boxSize_[1] << endl;
    return boxSize_;
  }

  void readPackingFromDirectory(string dirName, long numParticles_, long nDim_) {
    std::vector<double> pPos_(numParticles_ * nDim_);
    std::vector<double> pRad_(numParticles_);

    if(nDim_ == 2) {
      pPos_ = read2DFile(dirName + "pos.dat", numParticles_);
    } else if(nDim_ == 3) {
      pPos_ = read3DFile(dirName + "pos.dat", numParticles_);
    } else {
      cout << "FileIO::readParticlePackingFromDirectory: only dimensions 2 and 3 are allowed!" << endl;
    }
    sp_->setPositions(pPos_);
    pRad_ = read1DFile(dirName + "rad.dat", numParticles_);
    sp_->setRadii(pRad_);

    // set box dimensions
    std::vector<double> boxSize_(nDim_);
    boxSize_ = read1DFile(dirName + "boxSize.dat", nDim_);
    sp_->setBoxSize(boxSize_);
    boxSize_ = sp_->getBoxSize();
    if(nDim_ == 2) {
      cout << "FileIO::readPackingFromDirectory: phi: " << sp_->getPackingFraction() << " box-Lx: " << boxSize_[0] << ", Ly: " << boxSize_[1] << endl;
    } else if(nDim_ == 3) {
      cout << "FileIO::readPackingFromDirectory: phi: " << sp_->getPackingFraction() << " box-Lx: " << boxSize_[0] << ", Ly: " << boxSize_[1] << ", Lz: " << boxSize_[2] << endl;
    }
  }

  void readPBCPackingFromDirectory(string dirName, long numParticles_, long nDim_) {
    std::vector<double> pPos_(numParticles_ * nDim_);
    std::vector<double> pRad_(numParticles_);
    if(nDim_ == 2) {
      pPos_ = read2DFile(dirName + "pos.dat", numParticles_);
    } else if(nDim_ == 3) {
      pPos_ = read3DFile(dirName + "pos.dat", numParticles_);
    } else {
      cout << "FileIO::readPBCParticlePackingFromDirectory: only dimensions 2 and 3 are allowed!" << endl;
    }
    sp_->setPBCPositions(pPos_);
    pRad_ = read1DFile(dirName + "rad.dat", numParticles_);
    sp_->setRadii(pRad_);
    
    // set box dimensions
    std::vector<double> boxSize_(nDim_);
    boxSize_ = read1DFile(dirName + "boxSize.dat", nDim_);
    sp_->setBoxSize(boxSize_);
    boxSize_ = sp_->getBoxSize();
    if(nDim_ == 2) {
      cout << "FileIO::readPBCPackingFromDirectory: phi: " << sp_->getPackingFraction() << " box-Lx: " << boxSize_[0] << ", Ly: " << boxSize_[1] << endl;
    } else if(nDim_ == 3) {
      cout << "FileIO::readPBCPackingFromDirectory: phi: " << sp_->getPackingFraction() << " box-Lx: " << boxSize_[0] << ", Ly: " << boxSize_[1] << ", Lz: " << boxSize_[2] << endl;
    }
  }

  void saveParams(string dirName) {
    string fileParams = dirName + "params.dat";
    ofstream saveParams(fileParams.c_str());
    openOutputFile(fileParams);
    long numParticles = sp_->getNumParticles();
    saveParams << "numParticles" << "\t" << numParticles << endl;
    saveParams << "nDim" << "\t" << sp_->getNDim() << endl;
    saveParams << "sigma" << "\t" << 2 * sp_->getMeanParticleSize() << endl;
    saveParams << "epsilon" << "\t" << sp_->getEnergyCostant() << endl;
    saveParams << "dt" << "\t" << sp_->getTimeStep() << endl;
    saveParams << "phi" << "\t" << sp_->getPackingFraction() << endl;
    saveParams << "energy" << "\t" << sp_->getEnergy() / numParticles << endl;
    saveParams << "temperature" << "\t" << sp_->getTemperature() << endl;
    long numA = sp_->getTypeNumParticles();
    if(numA != numParticles) {
      saveParams << "numA" << "\t" << sp_->getTypeNumParticles() << endl;
    }
    saveParams.close();
  }

  void savePacking(string dirName) {
    saveParams(dirName);
    // save vectors
    long nDim = sp_->getNDim();
    save1DFile(dirName + "boxSize.dat", sp_->getBoxSize());
    save1DFile(dirName + "rad.dat", sp_->getRadii());
    save2DFile(dirName + "pos.dat", sp_->getPositions(), nDim);
    save2DFile(dirName + "vel.dat", sp_->getVelocities(), nDim);
    if(sp_->simControl.particleType == simControlStruct::particleEnum::active || sp_->simControl.particleType == simControlStruct::particleEnum::vicsek) {
      if(nDim == 2) {
        save1DFile(dirName + "angles.dat", sp_->getAngles());
      } else if(nDim == 3) {
        save2DFile(dirName + "angles.dat", sp_->getAngles(), nDim);
      } else {
        cout << "FileIO::savePacking: only dimensions 2 and 3 are allowed for angles!" << endl;
      }
    }
  }

  void savePBCPacking(string dirName) {
   saveParams(dirName);
    // save vectors
    long nDim = sp_->getNDim();
    save1DFile(dirName + "boxSize.dat", sp_->getBoxSize());
    save1DFile(dirName + "rad.dat", sp_->getRadii());
    save2DFile(dirName + "pos.dat", sp_->getPBCPositions(), nDim);
    save2DFile(dirName + "vel.dat", sp_->getVelocities(), nDim);
    if(sp_->simControl.particleType == simControlStruct::particleEnum::active || sp_->simControl.particleType == simControlStruct::particleEnum::vicsek) {
      if(nDim == 2) {
        save1DFile(dirName + "angles.dat", sp_->getAngles());
      } else if(nDim == 3) {
        save2DFile(dirName + "angles.dat", sp_->getAngles(), nDim);
      } else {
        cout << "FileIO::savePacking: only dimensions 2 and 3 are allowed for angles!" << endl;
      }
    }
  }

  void readVelocity(string dirName, long numParticles_, long nDim_) {
    std::vector<double> vel_(numParticles_ * nDim_);
    if(nDim_ == 2) {
      vel_ = read2DFile(dirName + "vel.dat", numParticles_);
    } else if(nDim_ == 3) {
      vel_ = read3DFile(dirName + "vel.dat", numParticles_);
    } else {
      cout << "FileIO::readVelocity: only dimensions 2 and 3 are allowed!" << endl;
    }
    sp_->setVelocities(vel_);
  }

  void readState(string dirName, long numParticles_, long nDim_, bool initAngles = false) {
    readVelocity(dirName, numParticles_, nDim_);
    if(initAngles == false) {
      if(sp_->simControl.particleType == simControlStruct::particleEnum::active || sp_->simControl.particleType == simControlStruct::particleEnum::vicsek) {
        std::vector<double> angle_(numParticles_);
        if(nDim_ == 2) {
          angle_ = read1DFile(dirName + "angles.dat", numParticles_);
        } else if(nDim_ == 3) {
          angle_.resize(numParticles_ * nDim_);
          angle_ = read3DFile(dirName + "angles.dat", numParticles_);
        } else {
          cout << "FileIO::readState: only dimensions 2 and 3 are allowed for angles!" << endl;
        }
        sp_->setAngles(angle_);
      }
    }
  }

  void saveState(string dirName) {
    save2DFile(dirName + "pos.dat", sp_->getPositions(), sp_->nDim);
    save2DFile(dirName + "vel.dat", sp_->getVelocities(), sp_->nDim);
    if(sp_->simControl.particleType == simControlStruct::particleEnum::active || sp_->simControl.particleType == simControlStruct::particleEnum::vicsek) {
      if(sp_->nDim == 2) {
        save1DFile(dirName + "angles.dat", sp_->getAngles());
      } else {
        save2DFile(dirName + "angles.dat", sp_->getAngles(), sp_->nDim);
      }
    }
  }

  void saveForces(string dirName) {
    save2DFile(dirName + "forces.dat", sp_->getForces(), sp_->nDim);
  }

  void saveNeighbors(string dirName) {
    if(sp_->simControl.neighborType == simControlStruct::neighborEnum::neighbor) {
      save2DIndexFile(dirName + "neighbors.dat", sp_->getNeighbors(), sp_->neighborListSize);
    }
  }

  void saveVicsekNeighbors(string dirName) {
    if(sp_->simControl.particleType == simControlStruct::particleEnum::vicsek) {
      save2DIndexFile(dirName + "vicsekNeighbors.dat", sp_->getVicsekNeighbors(), sp_->vicsekNeighborListSize);
    }
  }

  void saveDynamicsParams(string dirName, double damping) {
    string fileParams = dirName + "dynParams.dat";
    ofstream saveParams(fileParams.c_str());
    openOutputFile(fileParams);
    saveParams << "damping" << "\t" << damping << endl;
    if(sp_->simControl.particleType == simControlStruct::particleEnum::active) {
      double driving, taup;
      sp_->getSelfPropulsionParams(driving, taup);
      saveParams << "taup" << "\t" << taup << endl;
      saveParams << "f0" << "\t" << driving << endl;
    } else if(sp_->simControl.particleType == simControlStruct::particleEnum::vicsek) {
      double driving, taup, Jvicsek, Rvicsek;
      sp_->getVicsekParams(driving, taup, Jvicsek, Rvicsek);
      saveParams << "taup" << "\t" << taup << endl;
      saveParams << "f0" << "\t" << driving << endl;
      saveParams << "Rvicsek" << "\t" << Rvicsek << endl;
      saveParams << "Jvicsek" << "\t" << Jvicsek << endl;
    }
    saveParams.close();
  }

  void saveNoseHooverParams(string dirName) {
    double mass, damping;
    sp_->getNoseHooverParams(mass, damping);
    string fileParams = dirName + "nhParams.dat";
    ofstream saveParams(fileParams.c_str());
    openOutputFile(fileParams);
    saveParams << "mass" << "\t" << mass << endl;
    saveParams << "damping" << "\t" << damping << endl;
    if(sp_->simControl.particleType == simControlStruct::particleEnum::active) {
      double driving, taup;
      sp_->getSelfPropulsionParams(driving, taup);
      saveParams << "taup" << "\t" << taup << endl;
      saveParams << "f0" << "\t" << driving << endl;
    } else if(sp_->simControl.particleType == simControlStruct::particleEnum::vicsek) {
      double driving, taup, Jvicsek, Rvicsek;
      sp_->getVicsekParams(driving, taup, Jvicsek, Rvicsek);
      saveParams << "taup" << "\t" << taup << endl;
      saveParams << "f0" << "\t" << driving << endl;
      saveParams << "Rvicsek" << "\t" << Rvicsek << endl;
      saveParams << "Jvicsek" << "\t" << Jvicsek << endl;
    }
    saveParams.close();
  }

  void readNoseHooverParams(string dirName, double &mass, double &damping) {
    string fileParams = dirName + "nhParams.dat";
    ifstream readParams(fileParams.c_str());
    if (!readParams.is_open()) {
      cout << "Error: Unable to open file " << fileParams << " - setting default values" << endl;
      return;
    }
    string paramName;
    double paramValue;
    while (readParams >> paramName >> paramValue) {
      if(paramName == "mass") {
        mass = paramValue;
      } else if(paramName == "damping") {
        damping = paramValue;
      }
    }
    readParams.close();
    if(mass == 1 && damping == 1) {
      cout << "FileIO::saveNoseHooverParams: mass and damping are not saved in nhParams.dat! Setting mass and damping to 1" << endl;
    }
  }

  void saveDumpPacking(string dirName, long numParticles_, long nDim_, long timeStep_) {
    std::vector<double> pos(numParticles_ * nDim_), rad(numParticles_), boxSize(nDim_);
    pos = sp_->getPBCPositions();
    rad = sp_->getRadii();
    boxSize = sp_->getBoxSize();
    this->openOutputFile(dirName + "packing.dump");
    outputFile << "ITEM: TIMESTEP" << endl;
    outputFile << timeStep_ << endl;
    outputFile << "ITEM: NUMBER OF ATOMS" << endl;
    outputFile << numParticles_ << endl;
    if(nDim_ == 3) {
      outputFile << "ITEM: BOX BOUNDS pp pp fixed" << endl;
    } else {
      outputFile << "ITEM: BOX BOUNDS pp pp" << endl;
    }
    outputFile << 0 << "\t" << boxSize[0] << endl;
    outputFile << 0 << "\t" << boxSize[1] << endl;
    if(nDim_ == 3) {
      outputFile << 0 << "\t" << boxSize[2] << endl;
      outputFile << "ITEM: ATOMS id type radius xu yu zu" << endl;
    } else {
      outputFile << "ITEM: ATOMS id type radius xu yu" << endl;
    }
    int type = 1;
    for (long particleId = 0; particleId < numParticles_; particleId++) {
      if(particleId < sp_->numA) {
        type = 1;
      } else {
        type = 2;
      }
      //outputFile << particleId + 1 << "\t" << 1 << "\t" << particleId + 1 << "\t";
      //outputFile << rad[particleId] << "\t" << pos[particleId * nDim] << "\t" << pos[particleId * nDim + 1] << "\t" << pos[particleId * nDim + 2] << endl;
      outputFile << particleId + 1 << "\t" << type << "\t" << rad[particleId] << "\t";
      if(nDim_ == 3) {
        outputFile << pos[particleId * nDim_] << "\t" << pos[particleId * nDim_ + 1] << "\t" << pos[particleId * nDim_ + 2] << endl;
      } else {
        outputFile << pos[particleId * nDim_] << "\t" << pos[particleId * nDim_ + 1] << endl;
      }
    }
  }

};

#endif /* FILEIO_H */
