#define _USE_MATH_DEFINES

#include "Simulation.h"
#include <iostream>
#include <algorithm>
#include <omp.h>
#include <malloc.h>
#include <cmath>


using namespace std;

Simulation::Simulation(mt19937 genMain, uint nSimEntered)
{
    numSim = nSimEntered;
    gen = genMain;
}

Simulation::~Simulation() {}


