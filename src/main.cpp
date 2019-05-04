#include <iostream>
#include <vector>
#include <array>
#include <random>
#include <time.h>
#include "Simulation.h"

using uint = unsigned int;
using sint = short int;
using usint = unsigned short int;


//const usint nn = 7;
//const uint nn_arr[nn] = { 100, 250, 550, 1100, 3300, 5600, 12200 };
const usint nn = 1;
const uint nn_arr[nn] = { 500 };


int main()
{
	uint n_sim;
	std::cout << "Enter the number of simulations: ";
	std::cin >> n_sim;
	std::cout << std::endl;

	clock_t tStart = clock();

	std::mt19937 gen;
	for(usint nn_iter = 0; nn_iter < nn; nn_iter++)
	{
		for(uint sim_iter = 0; sim_iter < n_sim; sim_iter++)
		{
			gen.seed(time(0) + rand());
			Simulation Sim(gen, nn_arr[nn_iter]);
			Sim.f_nGen();
			Sim.f_GenerateXY();
			Sim.f_FillGraph();
            std::cout << "\tIteration:\t" << sim_iter + 1 << std::endl;
		}
	}

	std::cout << "Done! Execution time: " << (float)(clock() - tStart)/CLOCKS_PER_SEC << " sec." <<  std::endl;
	return 0;
}
