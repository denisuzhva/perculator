#define _USE_MATH_DEFINES
#include "Simulation.h"
#include <iostream>
#include <algorithm>
#include <omp.h>


Simulation::Simulation(std::mt19937 genMain, float eta_mean_main, uint n_sim_main) : gen(genMain), eta_mean(eta_mean_main), n_sim(n_sim_main)
{
    v_S.resize(n_sim);
    N_mean = (int)(eta_mean * S_0 / stringSigma);

    for(uint sim_iter = 0; sim_iter < n_sim; sim_iter++)
    {
        //tLog = clock();
        f_nGen();
        //std::cout << "\t f_nGen: \t" << (float)(clock() - tLog)/CLOCKS_PER_SEC << std::endl;

        //tLog = clock();
        f_GenerateXY();
        //std::cout << "\t f_GenerateXY: \t" << (float)(clock() - tLog)/CLOCKS_PER_SEC << std::endl;

        //tLog = clock();
        f_FillGraph();
        //std::cout << "\t f_FillGraph: \t" << (float)(clock() - tLog)/CLOCKS_PER_SEC << std::endl;

        //tLog = clock();
        //f_FindConnComp();
        //std::cout << "\t f_FindConnComp: \t" << (float)(clock() - tLog)/CLOCKS_PER_SEC << std::endl;

        //tLog = clock();
        //f_FindMulPtFB(sim_iter);
        //std::cout << "\t f_FindMulPtFB: \t" << (float)(clock() - tLog)/CLOCKS_PER_SEC << std::endl;

        //std::cout << "\nIteration:\t" << sim_iter + 1 << "\tcomplete" << std::endl;
        //std::cout << "*********************************\n";
    }
}


Simulation::~Simulation() { std::cout << "\nSim Ende\n"; }


//Math
inline usint Simulation::f_in_PDF_N(usint N_mean)
{
    //uniform_int_distribution<usint> disN(N_mean - (int)(pow((float)(N_mean), N_disp) + 0.5), N_mean + (int)(pow((float)(N_mean), N_disp) + 0.5));
    //uniform_int_distribution<usint> disN(N_mean - (int)((float)(N_mean)/3 + 0.5), N_mean + (int)((float)(N_mean)/3 + 0.5));
    //uniform_int_distribution<usint> disN(0, 2*N_mean);
    std::poisson_distribution<usint> disN(N_mean);
    return disN(gen);
    //return N_mean;
}


inline float Simulation::f_in_distXY(float xi, float xj, float yi, float yj, usint flag)
{
    return (sqrt((xi - xj) * (xi - xj) + (yi - yj) * (yi - yj)) - flag * rs);
}


inline float Simulation::f_in_PDF(float x, float y)
{
    return (2*(-x*x - y*y + 1)/3.14); // CHANGE ACCORDING TO R
}


// Initialize
void Simulation::f_nGen()
{
    N = f_in_PDF_N(N_mean);
    //N = N_mean;
    //std::cout << "Current string number:\t" << N << "\tMean:\t" << N_mean << std::endl;
}


void Simulation::f_GenerateXY()
{
    float ang, rad, xUni, yUni;
    std::uniform_real_distribution<float> disRad(0.0, 1.0);
    std::uniform_real_distribution<float> disAng(0.0, 2*M_PI);

    v_x.resize(N);
    v_y.resize(N);
    
    for(uint i = 0; i < N; i++)
    {
        rad = R*sqrt(disRad(gen));
        ang = disAng(gen);
        xUni = rad * cos(ang);
        yUni = rad * sin(ang);
        //cout << "X: " << xUni << " Y: " << yUni << endl;
        //cout << "Rad: " << sqrt(xUni*xUni + yUni*yUni) << " Ang: " << atan(yUni/xUni) * 180 / M_PI << " AngPrinc: " << ang * 180 / M_PI << endl;
        v_x[i] = xUni;
        v_y[i] = yUni;
    }
}


void Simulation::f_FillGraph()
{
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
void Simulation::f_FindConnComp() // used in constructor
{
    v_compData.clear(); // clear the vcm
    v_comp.clear(); // clear the temp cluster array
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


void Simulation::f_dfs(usint v)
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


/// cluster analysis
void Simulation::f_FindMulPtFB(uint sim_iter)
{
    float borderCoordAll[4]; // for the coordinates of the border of a cluster area: xL xR yD yU
    std::vector<float> v_xObs, v_yObs;
    v_xObs.clear();
    v_yObs.clear();
    float overallArea, MCDist, MCPointX, MCPointY, MCStep, n_k, S_k, eta_k, S_i;
    usint MCNThrown, MCNAll, N_k;

    S_i = 0;

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

            borderCoordAll[0] = *std::max_element(v_xObs.begin(), v_xObs.end()) + rs - 0.01; // x max [0]
            borderCoordAll[1] = *std::min_element(v_xObs.begin(), v_xObs.end()) - rs + 0.01; // x min [1]
            borderCoordAll[2] = *std::max_element(v_yObs.begin(), v_yObs.end()) + rs - 0.01; // y max [2]
            borderCoordAll[3] = *std::min_element(v_yObs.begin(), v_yObs.end()) - rs + 0.01; // y min [3]

            //MCPointX = borderCoordAll[1] + 0.003; // CHANGE HERE
            //MCPointY = borderCoordAll[3] + 0.003;
            MCNThrown = 0; // overall
            MCNAll = 0; // in the string
            MCDist = 0;
            MCStep = R*0.01*sqrt(sqrt(sqrt(N_k))); // triple
            //MCStep = R*0.01*sqrt(sqrt(sqrt(sqrt(N_k)))); // quadruple
            //MCStep = 0.01;

            //#pragma omp parallel for
            for(MCPointY = borderCoordAll[3]; MCPointY <= borderCoordAll[2]; MCPointY += MCStep)
            {
                for(MCPointX = borderCoordAll[1]; MCPointX <= borderCoordAll[0]; MCPointX += MCStep)
                {
                    MCNThrown++;
                    //#pragma omp parallel for
                    for(usint i = 0; i < N_k; i++)
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

            S_k = overallArea * (1 + 0.087 / N_k) * MCNAll / MCNThrown; // in order to approximate real string size at low N limit
            //S_k = overallArea * MCNAll / MCNThrown;

            S_i += S_k;
        }
    }
    
    v_S[sim_iter] = S_i;
}


// Out
void Simulation::f_sCalc()
{
	data_S.open("data_S.txt", std::ios_base::app);
    data_S << std::endl << N_mean << "\t\t" << eta_mean << "\t\t";

    S_Av = 0;
    S2_Av = 0;
    S_Disp = 0;
    S_Omega = 0;
    for(uint ittt = 0; ittt < v_S.size(); ittt++)
    {
        S_Av += (float) v_S[ittt];
    }
    S_Av /= n_sim;

    for(uint ittt = 0; ittt < v_S.size(); ittt++)
    {
        S2_Av += (float) v_S[ittt]*v_S[ittt];
    }
    S2_Av /= n_sim;

    S_Disp = S2_Av - S_Av*S_Av;
    S_Omega = S_Disp / S_Av;
    S_Omega /= stringSigma; // normed 

    data_S << S_Av << "\t" << S_Disp << "\t" << S_Omega << "\t";

    data_S.close();
}