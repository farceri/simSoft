# simSoft
Code for parallel MD simulations written in C++, OpenMP and Python.
Available integrators: NVE, Nose-Hoover and Langevin dynamics
Available potentials: harmonic, Lennard-Jones (LJ), Weeks-Chandler-Andersen (WCA), binary LJ, LJ-WCA, LJ+-

Usage
1. git clone
2. mkdir bin; make
3. change executable in makefile to compile different scripts (compress.cpp, test.cpp, runNVE.cpp)

Ex. 2D packing of polydisperse disks with packing fraction 0.72
![test](https://github.com/user-attachments/assets/8165be0f-c535-47a7-a537-fbd841e3a69a)
