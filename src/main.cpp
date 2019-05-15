#include <iostream>
#include "Simulation.h"
#include <time.h>

using namespace std;

int main()
{
    int numSim;
    cout << "Number of simulations (number of times the experiment will repeat). The more experiments, the more science you make: " << endl;
    cin >> numSim;
    //int numSim = 10;
    clock_t tStart = clock();
    for(int i = 0; i < numSim; i++)
        Simulation Sim(i);

    cout << "Done! Execution time: " << (double)(clock() - tStart)/CLOCKS_PER_SEC << " sec. Thanks for helping CERN, you've made a real physical research. Please, email the data_*.txt files at epirtochol@gmail.com" << endl;
    return 0;
}
