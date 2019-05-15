#define _USE_MATH_DEFINES
#include "Simulation.h"
#include <iostream>
#include <algorithm>
#include <omp.h>


Simulation::Simulation(std::mt19937 genMain, uint n_mean_main) : gen(genMain), N_mean(n_mean_main) {}


Simulation::~Simulation() {}


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
    return (2*(-x*x - y*y + 1)/3.14);
}


// Initialize
void Simulation::f_nGen()
{
    N = f_in_PDF_N(N_mean);
    //N = N_mean;
    std::cout << "Current string number:\t" << N << "\tMean:\t" << N_mean << std::endl;
}


void Simulation::f_GenerateXY()
{
    float ang, rad, xUni, yUni;
    std::uniform_real_distribution<float> disRad(0.0, 1.0);
    std::uniform_real_distribution<float> disAng(0.0, 2*M_PI);

    v_x = std::vector<float>(N);
    v_y = std::vector<float>(N);

    for(uint i = 0; i < N; i++)
    {
        rad = R * sqrt(disRad(gen));
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

    #pragma omp parallel for
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
    //v_compData.clear(); // clear the vcm
    //v_comp.clear(); // clear the temp cluster array
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


/// Cluster analysis
/*
void Simulation::f_FindMulPtFB()
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

            v_xObs = std::vector<float>(N_k);
            v_yObs = std::vector<float>(N_k);

            for(usint obsIter = 0; obsIter < v_compData[clusIter].size(); obsIter++)
            {
                v_xObs[obsIter] = v_x[v_compData[clusIter][obsIter]];
                v_yObs[obsIter] = v_y[v_compData[clusIter][obsIter]];
            }

            borderCoordAll[0] = *std::max_element(v_xObs.begin(), v_xObs.end()) + rs - 0.01; // x max [0]
            borderCoordAll[1] = *std::min_element(v_xObs.begin(), v_xObs.end()) - rs + 0.01; // x min [1]
            borderCoordAll[2] = *std::max_element(v_yObs.begin(), v_yObs.end()) + rs - 0.01; // y max [2]
            borderCoordAll[3] = *std::min_element(v_yObs.begin(), v_yObs.end()) - rs + 0.01; // y min [3]

            MCNThrown = 0; // overall
            MCNAll = 0; // in the string
            MCDist = 0;
            //MCStep = R*0.01*sqrt(sqrt(sqrt(N_k))); // triple
            //MCStep = R*0.01*sqrt(sqrt(sqrt(sqrt(N_k)))); // quadruple
            MCStep = R*0.01;

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

            //S_k = overallArea * (1 + 0.087 / N_k) * MCNAll / MCNThrown; // in order to approximate real string size
            S_k = overallArea * MCNAll / MCNThrown;

            n_k = sqrt(N_k * S_k / stringSigma);
            std::poisson_distribution<uint> disPois_nn(n_k);
            nF_k = disPois_nn(gen);
            nB_k = disPois_nn(gen);
            nF_i += (float)nF_k;
            nB_i += (float)nB_k;

            eta_k = N_k * stringSigma / S_k;
            sumPtF_av += (float)nF_k * sqrt(sqrt(eta_k));
            sumPtB_av += (float)nB_k * sqrt(sqrt(eta_k));
            sumPtF_disp += (float)nF_k * sqrt(eta_k);
            sumPtB_disp += (float)nB_k * sqrt(eta_k);
        }
    }
    pF_i_av = sumPtF_av / nF_i;
    pB_i_av = sumPtB_av / nB_i;
    pF_i_disp = sumPtF_disp / (nF_i * nF_i);
    pB_i_disp = sumPtB_disp / (nB_i * nB_i);
    std::normal_distribution<float> disNorm_ppF(pF_i_av, ptGamma*sqrt(pF_i_disp));
    std::normal_distribution<float> disNorm_ppB(pB_i_av, ptGamma*sqrt(pB_i_disp));
    pF_i = disNorm_ppF(gen);
    pB_i = disNorm_ppB(gen);
}
*/


void Simulation::f_FindMulPtFB()
{
    uint nF_k, nB_k;

    nF_i = 0;
    nB_i = 0;
    pF_i = 0;
    pB_i = 0;
    sumPtF_av = 0;
    sumPtB_av = 0;
    sumPtF_disp = 0;
    sumPtB_disp = 0;

    std::poisson_distribution<int> disPois(N);
    nF_i = disPois(gen);
    nB_i = disPois(gen);

    sumPtF_av += (float)nF_i;
    sumPtB_av += (float)nB_i;
    sumPtF_disp += (float)nF_i;
    sumPtB_disp += (float)nB_i;

    pF_i_av = 1;
    pB_i_av = 1;
    pF_i_disp = 1 / nF_i;
    pB_i_disp = 1 / nB_i;
    std::normal_distribution<float> disNorm_ppF(pF_i_av, ptGamma*sqrt(pF_i_disp));
    std::normal_distribution<float> disNorm_ppB(pB_i_av, ptGamma*sqrt(pB_i_disp));
    pF_i = disNorm_ppF(gen);
    pB_i = disNorm_ppB(gen);
}


// Out
float* Simulation::f_returnNP()
{
    static float np_arr[4];
    np_arr[0] = nF_i;
    np_arr[1] = nB_i;
    np_arr[2] = pF_i;
    np_arr[3] = pB_i;
    //std::cout << "\nTrue nF:\t" << nF_i << "\t";
    //std::cout << "\nTrue nB:\t" << nB_i << "\t";
    //std::cout << "\nTrue pF:\t" << pF_i << "\t";
    //std::cout << "\nTrue pB:\t" << pB_i << "\t";
    return np_arr;
}