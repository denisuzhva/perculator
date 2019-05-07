#include <iostream>
#include <array>
#include <random>
#include "Simulation.h"


using uint = unsigned int;
using sint = short int;
using usint = unsigned short int;


int main()
{
	uint n_sim;
	std::cout << "Enter the number of simulations: ";
	std::cin >> n_sim;

	clock_t tStart = clock();

	
	const usint nn = 4;
	const float eta_arr[nn] = {0.05, 0.1, 0.15, 0.2};
	//const float eta_arr[nn] = {0.2};
	//const uint nn_arr[nn] = {6600};
	//const uint nn_arr[nn] = {100, 250, 550, 1100, 3300, 5600, 12200};
	//uint nn_arr[nn];
	std::mt19937 gen;

	for(usint nn_iter = 0; nn_iter < nn; nn_iter++)
	{
		std::cout << "\nN_mean:\t" << eta_arr[nn_iter] << std::endl;

		gen.seed(time(0) + rand());
		Simulation Sim(gen, eta_arr[nn_iter], n_sim);	
		Sim.f_sCalc();
	}


	std::cout << "\nDone! Execution time: " << (float)(clock() - tStart)/CLOCKS_PER_SEC << " sec." <<  std::endl;
	getchar();
	return 0;
}
