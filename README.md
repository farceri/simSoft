# simSoft
Code for parallel MD simulations written in C++, OpenMP and Python.
Available integrators: NVE, Nose-Hoover and Langevin dynamics
Available potentials: harmonic, Lennard-Jones (LJ), Weeks-Chandler-Andersen (WCA), binary LJ, LJ-WCA, LJ+-

Usage
1. git clone
2. mkdir bin; make
3. change executable in makefile to compile different scripts (compress.cpp, test.cpp, runNVE.cpp)

Ex. 2D packing of polydisperse disks with packing fraction 0.72
![test](https://github.com/user-attachments/assets/d25c6193-3eeb-4871-95bb-578d2aa3a6ec)
