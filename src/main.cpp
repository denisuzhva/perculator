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

	
	//const usint nn = 5;
	//const uint nn_arr[nn] = {100, 250, 550, 1100, 3300};
	const usint nn = 1;
	//const uint nn_arr[nn] = {5600};
	const float eta_arr[nn] = {5};
	std::mt19937 gen;
	std::vector<uint> nF_vec, nB_vec;
	std::vector<float> pF_vec, pB_vec;
	nF_vec.resize(n_sim);
	nB_vec.resize(n_sim);
	pF_vec.resize(n_sim);
	pB_vec.resize(n_sim);
	float nFnB_av, nF_av, nB_av, nF2_av, nF_av2; // averaging for the multiplicities
    float pFpB_av, pF_av, pB_av, pF2_av, pF_av2; // averaging for the pt
    float b_nn, b_pp; //


	for(usint nn_iter = 0; nn_iter < nn; nn_iter++)
	{

		std::cout << "\neta_mean:\t" << eta_arr[nn_iter] << std::endl;

		nFnB_av = 0;
		nF_av = 0;
		nB_av = 0;
		nF2_av = 0;
		nF_av2 = 0;
		pFpB_av = 0;
		pF_av = 0;
		pB_av = 0;
		pF2_av = 0;
		pF_av2 = 0;

		for(uint sim_iter = 0; sim_iter < n_sim; sim_iter++)
		{
			gen.seed(time(0) + rand());
			Simulation Sim(gen, eta_arr[nn_iter]);
			Sim.f_nGen();
			Sim.f_GenerateXY();
			Sim.f_FillGraph();
			Sim.f_FindConnComp();

			//tLog = clock();
			Sim.f_FindMulPtFB();
			//std::cout << "\t f_FindMulPtFB: \t" << (float)(clock() - tLog)/CLOCKS_PER_SEC << std::endl;


			//std::cout << "\nTrue nF:\t" << nF_vec[sim_iter] << "\t";
    		//std::cout << "\nTrue nB:\t" << nB_vec[sim_iter] << "\t";

            //std::cout << "\nIteration:\t" << sim_iter + 1 << "\tcompleted" << std::endl;
			//std::cout << "*********************************\n";
		}		
	}
	

	std::cout << "\nDone! Execution time: " << (float)(clock() - tStart)/CLOCKS_PER_SEC << " sec." <<  std::endl;
	getchar();
	return 0;
}