#include <iostream>
#include <vector>
#include <array>
#include <random>
#include <time.h>
#include <fstream>
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
	clock_t tLog;

	
	const usint nn = 7;
	const uint nn_arr[nn] = {100, 250, 550, 1100, 3300, 5600, 12200};
	//const usint nn = 1;
	//const uint nn_arr[nn] = {5600};
	std::mt19937 gen;
	std::vector<float> pF_vec, pB_vec;
	pF_vec.resize(n_sim);
	pB_vec.resize(n_sim);
	std::ofstream data_pF_i, data_pB_i;
	float *np;

	data_pF_i.open("data_pF_i.txt", std::ios_base::app);
	data_pB_i.open("data_pB_i.txt", std::ios_base::app);


	for(usint nn_iter = 0; nn_iter < nn; nn_iter++)
	{
		std::cout << "\nN_mean:\t" << nn_arr[nn_iter] << std::endl;
		data_pF_i << std::endl << nn_arr[nn_iter] << "\t\t";
		data_pB_i << std::endl << nn_arr[nn_iter] << "\t\t";

		for(uint sim_iter = 0; sim_iter < n_sim; sim_iter++)
		{
			gen.seed(time(0) + rand());
			Simulation Sim(gen, nn_arr[nn_iter]);
			Sim.f_nGen();
			Sim.f_GenerateXY();
			Sim.f_FillGraph();
			Sim.f_FindConnComp();

			//tLog = clock();
			Sim.f_FindMulPtFB();
			//std::cout << "\t f_FindMulPtFB: \t" << (float)(clock() - tLog)/CLOCKS_PER_SEC << std::endl;

			np = Sim.f_returnNP();
			pF_vec[sim_iter] = np[0];
			pB_vec[sim_iter] = np[1];

			//std::cout << "\nTrue nF:\t" << nF_vec[sim_iter] << "\t";
    		//std::cout << "\nTrue nB:\t" << nB_vec[sim_iter] << "\t";

    		data_pF_i << pF_vec[sim_iter] << "\t";
    		data_pB_i << pB_vec[sim_iter] << "\t";

            std::cout << "\nIteration:\t" << sim_iter + 1 << "\tcompleted" << std::endl;
			std::cout << "*********************************\n";
		}
	}
	

    data_pF_i.close();
    data_pB_i.close();

	std::cout << "\nDone! Execution time: " << (float)(clock() - tStart)/CLOCKS_PER_SEC << " sec." <<  std::endl;
	getchar();
	return 0;
}