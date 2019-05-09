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

	
	//const usint nn = 7;
	//const uint nn_arr[nn] = {100, 250, 550, 1100, 3300, 5600, 12200};
	const usint nn = 1;
	const uint nn_arr[nn] = {5600};
	std::mt19937 gen;
	std::vector<uint> nF_vec, nB_vec;
	std::vector<float> pF_vec, pB_vec;
	nF_vec.resize(n_sim);
	nB_vec.resize(n_sim);
	pF_vec.resize(n_sim);
	pB_vec.resize(n_sim);
	std::ofstream data_nF_i, data_nB_i, data_pF_i, data_pB_i, data_b;
	float nFnB_av, nF_av, nB_av, nF2_av, nF_av2; // averaging for the multiplicities
    float pFpB_av, pF_av, pB_av, pF2_av, pF_av2; // averaging for the pt
    float b_nn, b_pp; //
	float *np;

	data_nF_i.open("data_nF_i.txt", std::ios_base::app);
	data_nB_i.open("data_nB_i.txt", std::ios_base::app);
	data_pF_i.open("data_pF_i.txt", std::ios_base::app);
	data_pB_i.open("data_pB_i.txt", std::ios_base::app);
	data_b.open("data_b.txt", std::ios_base::app);


	for(usint nn_iter = 0; nn_iter < nn; nn_iter++)
	{
		std::cout << "\nN_mean:\t" << nn_arr[nn_iter] << std::endl;
		data_nF_i << std::endl << nn_arr[nn_iter] << "\t\t";
		data_nB_i << std::endl << nn_arr[nn_iter] << "\t\t";
		data_pF_i << std::endl << nn_arr[nn_iter] << "\t\t";
		data_pB_i << std::endl << nn_arr[nn_iter] << "\t\t";
		data_b << std::endl << nn_arr[nn_iter] << "\t\t";

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
			nF_vec[sim_iter] = np[0];
			nB_vec[sim_iter] = np[1];
			pF_vec[sim_iter] = np[2];
			pB_vec[sim_iter] = np[3];

			std::cout << "\nTrue nF:\t" << nF_vec[sim_iter] << "\t";
    		std::cout << "\nTrue nB:\t" << nB_vec[sim_iter] << "\t";

			data_nF_i << nF_vec[sim_iter] << "\t";
    		data_nB_i << nB_vec[sim_iter] << "\t";
    		data_pF_i << pF_vec[sim_iter] << "\t";
    		data_pB_i << pB_vec[sim_iter] << "\t";

            std::cout << "\nIteration:\t" << sim_iter + 1 << "\tcompleted" << std::endl;
			std::cout << "*********************************\n";
		}

		// nn
		for(uint i = 0; i < n_sim; i++)
			nFnB_av += (float)(nF_vec[i]*nB_vec[i]);
		nFnB_av /= n_sim;

		for(uint i = 0; i < n_sim; i++)
			nF_av += (float)(nF_vec[i]);
		nF_av /= n_sim;

		for(uint i = 0; i < n_sim; i++)
			nB_av += (float)(nB_vec[i]);
		nB_av /= n_sim;

		for(uint i = 0; i < n_sim; i++)
			nF2_av += (float)(nF_vec[i]*nF_vec[i]);
		nF2_av /= n_sim;

		nF_av2 = (float)(nF_av*nF_av);

		// pp
		for(uint i = 0; i < n_sim; i++)
			pFpB_av += (float)(pF_vec[i]*pB_vec[i]);
		pFpB_av /= n_sim;

		for(uint i = 0; i < n_sim; i++)
			pF_av += (float)(pF_vec[i]);
		pF_av /= n_sim;

		for(uint i = 0; i < n_sim; i++)
			pB_av += (float)(pB_vec[i]);
		pB_av /= n_sim;

		for(uint i = 0; i < n_sim; i++)
			pF2_av += (float)(pF_vec[i]*pF_vec[i]);
		pF2_av /= n_sim;

		pF_av2 = (float)(pF_av*pF_av);

		// b
		b_nn = (nFnB_av - nF_av*nB_av) / (nF2_av - nF_av2);
    	b_pp = (pFpB_av - pF_av*pB_av) / (pF2_av - pF_av2);

    	data_b << b_nn << "\t" << b_pp << "\t";
	}
	

	data_nF_i.close();
    data_nB_i.close();
    data_pF_i.close();
    data_pB_i.close();
    data_b.close();

	std::cout << "\nDone! Execution time: " << (float)(clock() - tStart)/CLOCKS_PER_SEC << " sec." <<  std::endl;
	getchar();
	return 0;
}