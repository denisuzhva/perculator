#include <iostream>
#include "Simulation.h"
#include <random>
#include <time.h>

using uint = unsigned int;
using sint = short int;
using usint = unsigned short int;



int main()
{
	const uint nn = 7;
	uint nn_arr[nn] = {100, 250, 550, 1100, 3300, 5600, 12200};

	uint n_sim;
	std::cout << "Enter the number of simulations: ";
	std::cin >> n_sim;
	std::cout << std::endl;
	

	clock_t tStart = clock();

	mt19937* gen(time(0));
	for(uint nn_iter = 0; nn_iter < nn; nn_iter++)
	{
		for(uint sim_iter = 0; sim_iter < n_sim; sim_iter++)
		{
			
		}
	}

	std::cout << "Done! Execution time: " << (float)(clock() - tStart)/CLOCKS_PER_SEC << " sec." <<  std::endl;
	return 0;
}
