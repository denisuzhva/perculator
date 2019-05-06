#include <iostream>
#include <vector>
#include <algorithm>
#include <random>
#include <time.h>
#include <fstream>

#define _USE_MATH_DEFINES
#define M_PI       3.14159265358979323846 


using uint = unsigned int;
using sint = short int;
using usint = unsigned short int;


//const usint nn = 7;
//const usint nn_arr[nn] = {100, 250, 550, 1100, 3300, 5600, 12200};
const usint nn = 1;
const uint nn_arr[nn] = {15000};

constexpr static float rs = 0.03, stringSigma = 3.14*rs*rs; // string's radius (def. 0.03) and area
constexpr static float ptGammaSquared = 0.5; // proportionality coefficient in pt distribution

std::mt19937 gen(time(0));
usint N = 0;
std::vector<float> v_x; // vector of x coordinates
std::vector<float> v_y; // vector of y coordinate
std::vector<std::vector<usint>> g_connGraph; // graph for conn comps
std::vector<bool> v_used; // vector of used comps in dfs
std::vector<usint> v_comp; // vector for temp conn comps
std::vector<std::vector<usint>> v_compData; // matrix for conn comps
float nF_i, nB_i, pF_i, pB_i; // F, B for 1 simulation
float sumPtF_av, sumPtB_av, sumPtF_disp, sumPtB_disp; // sum in av & disp pt calc
float pF_i_av, pB_i_av, pF_i_disp, pB_i_disp; // av and disp pt calc


//Math
inline usint f_in_PDF_N(usint N_mean)
{
    //uniform_int_distribution<usint> disN(N_mean - (int)(pow((float)(N_mean), N_disp) + 0.5), N_mean + (int)(pow((float)(N_mean), N_disp) + 0.5));
    //uniform_int_distribution<usint> disN(N_mean - (int)((float)(N_mean)/3 + 0.5), N_mean + (int)((float)(N_mean)/3 + 0.5));
    //uniform_int_distribution<usint> disN(0, 2*N_mean);
    std::poisson_distribution<usint> disN(N_mean);
    return disN(gen);
    //return N_mean;
}


inline float f_in_distXY(float xi, float xj, float yi, float yj, usint flag)
{
    return (sqrt((xi - xj) * (xi - xj) + (yi - yj) * (yi - yj)) - flag * rs);
}


inline float f_in_PDF(float x, float y)
{
    return (2*(-x*x - y*y + 1)/3.14);
}


// Initialize
void f_nGen(usint N_mean)
{
	N = f_in_PDF_N(N_mean);
    //std::cout << "Current string number:\t" << N << "\tMean:\t" << N_mean << std::endl;
}


void f_GenerateXY()
{
    float ang, rad, xUni, yUni;
    std::uniform_real_distribution<float> disRad(0.0, 1.0);
    std::uniform_real_distribution<float> disAng(0.0, 2*M_PI);

	v_x.clear();
	v_y.clear();

    v_x.resize(N);
    v_y.resize(N);

    for(uint i = 0; i < N; i++)
    {
        rad = sqrt(disRad(gen));
        ang = disAng(gen);
        xUni = rad * cos(ang);
        yUni = rad * sin(ang);
        //cout << "X: " << xUni << " Y: " << yUni << endl;
        //cout << "Rad: " << sqrt(xUni*xUni + yUni*yUni) << " Ang: " << atan(yUni/xUni) * 180 / M_PI << " AngPrinc: " << ang * 180 / M_PI << endl;
        v_x[i] = xUni;
        v_y[i] = yUni;
    }

}


void f_FillGraph()
{
	g_connGraph.clear();
    g_connGraph.resize(N, std::vector<usint>(N));

    for(usint i = 0; i < N; i++)
        for(usint j = 0; j < N; j++)
        {
            g_connGraph[i][j] = 65535; // unsigned short int is used
            //cout << -1;
        }

    //#pragma omp parallel for
    for(usint i = 0; i < N; i++)
    {
        for(usint j = 0; j < N; j++)
        {
            if(i != j)
                if(f_in_distXY(v_x[i], v_x[j], v_y[i], v_y[j], 2) <= 0)
                {
                    g_connGraph[i][j] = j; // the graph made out of 16bits and numbers of all the centers of strings overlapped with j'th
                }
        }
    }
}


/// Conn comp calculating
void f_dfs(usint v)
{
    v_used[v] = true; // set v'th coordinate as being searched for near strings, overlapping with v'th
	v_comp.push_back(v); // push v'th string to a component being made
	for(usint i = 0; i < N; ++i)
	{
	    if(g_connGraph[v][i] < 65535) // do if not equal to -1 (if GCG is usint then < 65535)
        {
            usint to = g_connGraph[v][i]; // get a neighbor from the graph
            if(! v_used[to]) // check if a neighbor was not checked yet
                f_dfs(to); // recursively repeat dfs with all the found neighbors not checked yet (unused)
        }
	}
}


void f_FindConnComp() // used in constructor
{
    v_compData.clear(); // clear the vcm
    v_comp.clear(); // clear the temp cluster array
	v_used.clear();
    v_used.resize(N);
    for(usint i = 0; i < N; ++i)
		v_used[i] = false; // make all the coordinates unused
	for(usint i = 0; i < N; ++i)
		if(!v_used[i]) // if unused...
        {
            v_comp.clear();
			f_dfs(i); // go to the depth-first search
			if(!v_comp.empty())
            {
                v_compData.resize(i+1);
                v_compData[i] = v_comp; // write the calculated component into the vcm
            }
		}
}


/// cluster analysis
void f_FindMulPtFB()
{
    float borderCoordAll[4]; // for the coordinates of the border of a cluster area: xL xR yD yU
    std::vector<float> v_xObs, v_yObs;
    float overallArea, MCDist, MCPointX, MCPointY, MCStep, n_k, S_k, eta_k;
    usint MCNThrown, MCNAll, N_k;
    uint nF_k, nB_k;

    //cout << v_compData.size() << endl;

    nF_i = 0;
    nB_i = 0;
    pF_i = 0;
    pB_i = 0;
    sumPtF_av = 0;
    sumPtB_av = 0;
    sumPtF_disp = 0;
    sumPtB_disp = 0;

    for(usint clusIter = 0; clusIter < v_compData.size(); clusIter++)
    {
        if(v_compData[clusIter].size() != 0)
        {
            N_k = v_compData[clusIter].size();

            v_xObs.resize(N_k);
            v_yObs.resize(N_k);

            for(usint obsIter = 0; obsIter < v_compData[clusIter].size(); obsIter++)
            {
                v_xObs[obsIter] = v_x[v_compData[clusIter][obsIter]];
                v_yObs[obsIter] = v_y[v_compData[clusIter][obsIter]];
            }

            borderCoordAll[0] = *std::max_element(v_xObs.begin(), v_xObs.end()) + rs - 0.001; // x max [0]
            borderCoordAll[1] = *std::min_element(v_xObs.begin(), v_xObs.end()) - rs + 0.001; // x min [1]
            borderCoordAll[2] = *std::max_element(v_yObs.begin(), v_yObs.end()) + rs - 0.001; // y max [2]
            borderCoordAll[3] = *std::min_element(v_yObs.begin(), v_yObs.end()) - rs + 0.001; // y min [3]

            MCPointX = borderCoordAll[1] + 0.003;
            MCPointY = borderCoordAll[3] + 0.003;
            MCNThrown = 0; // overall
            MCNAll = 0; // in the string
            MCDist = 0;
            //MCStep = 0.01*sqrt(sqrt(sqrt(N_k))); // triple
            MCStep = 0.01*sqrt(sqrt(sqrt(sqrt(N_k)))); // quadruple
            //MCStep = 0.01;

            //#pragma omp parallel for
            for(MCPointY = borderCoordAll[3]; MCPointY <= borderCoordAll[2]; MCPointY += MCStep)
            {
                for(MCPointX = borderCoordAll[1]; MCPointX <= borderCoordAll[0]; MCPointX += MCStep)
                {
                    MCNThrown++;
                    //#pragma omp parallel for
                    for(usint i = 0; i < v_xObs.size(); i++) // no matter v_xObs or v_yObs size
                    {
                        MCDist = f_in_distXY(MCPointX, v_xObs[i], MCPointY, v_yObs[i], 0);
                        if(MCDist <= rs)
                        {
                            MCNAll++;
                            break;
                        }
                    }
                }
            }

            overallArea = (borderCoordAll[0] - borderCoordAll[1]) * (borderCoordAll[2] - borderCoordAll[3]);

            S_k = overallArea * (1 + 0.087 / N_k) * MCNAll / MCNThrown; // in order to approximate real string size
            //S_k = overallArea * MCNAll / MCNThrown;

            n_k = sqrt(N_k * S_k / stringSigma);
            std::poisson_distribution<uint> disPois_nn(n_k);
            nF_k = disPois_nn(gen);
            nB_k = disPois_nn(gen);
            nF_i += (float)nF_k;
            nB_i += (float)nB_k;

            eta_k = N_k * stringSigma / S_k;
            sumPtF_av += nF_k * sqrt(sqrt(eta_k));
            sumPtB_av += nF_k * sqrt(sqrt(eta_k));
            sumPtF_disp += nF_k * sqrt(eta_k);
            sumPtB_disp += nB_k * sqrt(eta_k);
        }
    }
    pF_i_av = sumPtF_av / nF_i;
    pB_i_av = sumPtB_av / nB_i;
    pF_i_disp = sumPtF_disp / (nF_i * nF_i);
    pB_i_disp = sumPtB_disp / (nB_i * nB_i);
    std::normal_distribution<float> disNorm_ptptF(pF_i_av, ptGammaSquared*sqrt(pF_i_disp));
    std::normal_distribution<float> disNorm_ptptB(pB_i_av, ptGammaSquared*sqrt(pB_i_disp));
    pF_i = disNorm_ptptF(gen);
    pB_i = disNorm_ptptB(gen);
}


int main()
{
	uint n_sim;
	std::cout << "Enter the number of simulations: ";
	std::cin >> n_sim;

	clock_t tStart = clock();

	std::vector<uint> nF_vec, nB_vec;
	std::vector<float> pF_vec, pB_vec;
	nF_vec.resize(n_sim);
	nB_vec.resize(n_sim);
	pF_vec.resize(n_sim);
	pB_vec.resize(n_sim);
	std::ofstream data_nF_i, data_nB_i, data_pF_i, data_pB_i, data_b;
	float nFnB_av, nF_av, nB_av, nF2_av, nF_av2; // averaging for the multiplicities
    float pFpB_av, pF_av, pB_av, pF2_av, pF_av2; // averaging for the pt
    float b_nn, b_pp;

	//data_nF_i.open("data_nF_i.txt", std::ios_base::app);
	//data_nB_i.open("data_nB_i.txt", std::ios_base::app);
	//data_pF_i.open("data_pF_i.txt", std::ios_base::app);
	//data_pB_i.open("data_pB_i.txt", std::ios_base::app);
	data_b.open("data_b.txt", std::ios_base::app);


	for(usint nn_iter = 0; nn_iter < nn; nn_iter++)
	{
		//std::cout << "\nN_mean:\t" << nn_arr[nn_iter];
		//data_nF_i << std::endl << nn_arr[nn_iter] << "\t\t";
		//data_nB_i << std::endl << nn_arr[nn_iter] << "\t\t";
		//data_pF_i << std::endl << nn_arr[nn_iter] << "\t\t";
		//data_pB_i << std::endl << nn_arr[nn_iter] << "\t\t";
		data_b << std::endl << nn_arr[nn_iter] << "\t\t";

		for(uint sim_iter = 0; sim_iter < n_sim; sim_iter++)
		{
			//gen.seed(time(0) + rand());
			f_nGen(nn_arr[nn_iter]);
			f_GenerateXY();
			f_FillGraph();
			f_FindConnComp();
			f_FindMulPtFB();
			nF_vec[sim_iter] = nF_i;
			nB_vec[sim_iter] = nB_i;
			pF_vec[sim_iter] = pF_i;
			pB_vec[sim_iter] = pB_i;

			//data_nF_i << nF_vec[sim_iter] << "\t";
    		//data_nB_i << nB_vec[sim_iter] << "\t";
    		//data_pF_i << pF_vec[sim_iter] << "\t";
    		//data_pB_i << pB_vec[sim_iter] << "\t";

            //std::cout << "\nIteration:\t" << sim_iter + 1 << "\tcomplete" << std::endl;
			//std::cout << "*********************************\n";
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
	

	//data_nF_i.close();
    //data_nB_i.close();
    //data_pF_i.close();
    //data_pB_i.close();
    data_b.close();

	std::cout << "\nDone! Execution time: " << (float)(clock() - tStart)/CLOCKS_PER_SEC << " sec." <<  std::endl;
	getchar();
	return 0;
}
