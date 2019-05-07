#include <iostream>
#include "include/Simulation.h"
#include <random>
#include <time.h>

using namespace std;

int main()
{
    uint nSim;
    //uint dispersionz[9] = {1, 3, 5, 10, 15, 20, 25, 50, 100};
    //float dispersionz[9] = {1.0, 0.9, 0.8, 0.7, 0.6, 0.5, 0.4, 0.3, 0.2};
    //float disp;
    cout << "Enter the number of simulations: ";
    cin >> nSim;
    cout << endl;
    clock_t tStart = clock();

    mt19937 gen(time(0));
    Simulation Sim(gen, nSim);

    cout << "Done! Execution time: " << (float)(clock() - tStart)/CLOCKS_PER_SEC << " sec." << endl;
    return 0;
}
