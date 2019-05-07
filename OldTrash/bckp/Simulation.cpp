#define _USE_MATH_DEFINES

#include "Simulation.h"
#include <iostream>
#include <algorithm>
#include <omp.h>
#include <malloc.h>
#include <cmath>
#include <time.h>


using namespace std;

Simulation::Simulation(mt19937 genMain, uint nSimEntered)
{
    numSim = nSimEntered;
    gen = genMain;
    for(usint iN = 0; iN < a_nRange.size(); iN++)
    {
        N_mean = a_nRange[iN];
        cout << "\nN_mean:\t" << N_mean << endl;

        nFnB_Av = 0;
        nF_Av = 0;
        nB_Av = 0;
        nF2_Av = 0;
        nF_Av2 = 0;
        v_nF_i.clear();
        v_nB_i.clear();

        pFpB_Av = 0;
        pF_Av = 0;
        pB_Av = 0;
        pF2_Av = 0;
        pF_Av2 = 0;
        v_pF_i.clear();
        v_pB_i.clear();

        data_nF_i.open("data_nF_i.txt", ios_base::app);
        data_nF_i << endl << N_mean << "\t\t";

        data_nB_i.open("data_nB_i.txt", ios_base::app);
        data_nB_i << endl << N_mean << "\t\t";

        data_pF_i.open("data_pF_i.txt", ios_base::app);
        data_pF_i << endl << N_mean << "\t\t";

        data_pB_i.open("data_pB_i.txt", ios_base::app);
        data_pB_i << endl << N_mean << "\t\t";

        data_b.open("data_b.txt", ios_base::app);
        data_b << endl << N_mean << "\t\t";

        for(uint iSim = 0; iSim < numSim; iSim++)
        {
            f_nGen();
            //cout << "\tIteration:\t" << iSim+1 << endl;
            f_GenerateXY();
            //f_GenerateXY_Neym();
            f_FillGraph();
            f_FindConnComp();

			tLog = clock();
            f_FindMulPtFB();
			std::cout << "\t f_FindMulPtFB: \t" << (float)(clock() - tLog)/CLOCKS_PER_SEC << std::endl;

            f_WriteData();
        }

        f_bCalc();

        data_nF_i.close();
        data_nB_i.close();
        data_pF_i.close();
        data_pB_i.close();
        data_b.close();
    }

}

Simulation::~Simulation() {}


/// initialization
void Simulation::f_nGen()
{
    //N = f_in_PDF_N(N_mean);
    N = N_mean;
    cout << "Current string number:\t" << N << "\tMean:\t" << N_mean << endl;
}


void Simulation::f_GenerateXY()
{
    v_x.clear();
    v_y.clear();
    float ang, rad, xUni, yUni;
    uniform_real_distribution<float> disRad(0.0, 1.0);
    uniform_real_distribution<float> disAng(0.0, 2*M_PI);

    for(uint i = 0; i < N; i++)
    {
        rad = sqrt(disRad(gen));
        ang = disAng(gen);
        xUni = rad * cos(ang);
        yUni = rad * sin(ang);
        //cout << "X: " << xUni << " Y: " << yUni << endl;
        //cout << "Rad: " << sqrt(xUni*xUni + yUni*yUni) << " Ang: " << atan(yUni/xUni) * 180 / M_PI << " AngPrinc: " << ang * 180 / M_PI << endl;
        v_x.push_back(xUni);
        v_y.push_back(yUni);
    }

}


void Simulation::f_GenerateXY_Neym()
{
    v_x.clear();
    v_y.clear();
    float xNeym, yNeym, z;
    usint i = 0;
    uniform_real_distribution<float> dis1(-1.0, 1.0);
    uniform_real_distribution<float> dis2(0, 1.0);
    while(i < N)
    {
        xNeym = dis1(gen);
        yNeym = dis1(gen);
        z = dis2(gen);
        if(z < f_in_PDF(xNeym, yNeym))
        {
            v_x.push_back(xNeym);
            v_y.push_back(yNeym);
            ++i;
        }
    }
}

void Simulation::f_FillGraph()
{
    g_connGraph.clear();

    g_connGraph.resize(N);
    for(usint i = 0; i < N; i++)
        g_connGraph[i].resize(N);


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
                    g_connGraph[i][j] = j; // the graph made out of NULL's and numbers of all the centers of strings overlapped with j'th

                }
        }
    }
}


/// conn comp calculating
void Simulation::f_FindConnComp() // used in constructor
{
    v_compData.clear(); // clear the vcm
    v_comp.clear(); // clear the temp cluster array
    v_used.clear(); // clear the trigger marks
    for(usint i = 0; i < N; ++i)
		v_used.push_back(false); // make all the coordinates unused
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


/// math
float Simulation::f_in_PDF(float x, float y)
{
    return (2*(-x*x - y*y + 1)/3.14);
}

usint Simulation::f_in_PDF_N(usint N_mean)
{
    //uniform_int_distribution<usint> disN(N_mean - (int)(pow((float)(N_mean), N_disp) + 0.5), N_mean + (int)(pow((float)(N_mean), N_disp) + 0.5));
    //uniform_int_distribution<usint> disN(N_mean - (int)((float)(N_mean)/3 + 0.5), N_mean + (int)((float)(N_mean)/3 + 0.5));
    //uniform_int_distribution<usint> disN(0, 2*N_mean);
    poisson_distribution<usint> disN(N_mean);
    return disN(gen);
    //return N_mean;
}

/*
double Simulation::f_in_PDF_mul(double mulPDF)
{
    return mulPDF;
}
*/

float Simulation::f_in_distXY(float xi, float xj, float yi, float yj, usint factor)
{
    return (sqrt((xi - xj) * (xi - xj) + (yi - yj) * (yi - yj)) - factor * rs);
}



/// cluster analysis
void Simulation::f_FindMulPtFB()
{
    float *borderCoordAll; // for the coordinates of the border of a cluster area: xL xR yD yU
    vector<float> v_xObs, v_yObs;
    v_xObs.clear();
    v_yObs.clear();
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

            v_xObs.clear();
            v_yObs.clear();
            for(usint obsIter = 0; obsIter < v_compData[clusIter].size(); obsIter++)
            {
                v_xObs.push_back(v_x[v_compData[clusIter][obsIter]]);
                v_yObs.push_back(v_y[v_compData[clusIter][obsIter]]);
            }

            borderCoordAll = new float [4];

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

            overallArea = (borderCoordAll[0] - borderCoordAll[1])*(borderCoordAll[2] - borderCoordAll[3]);

            S_k = overallArea * (1 + 0.087 / N_k) * MCNAll / MCNThrown; // in order to approximate real string size
            //S_k = overallArea * MCNAll / MCNThrown;

            n_k = sqrt(N_k * S_k / stringSigma);
            poisson_distribution<uint> disPois_nn(n_k);
            nF_k = disPois_nn(gen);
            nB_k = disPois_nn(gen);
            nF_i += nF_k;
            nB_i += nB_k;

            eta_k = N_k * stringSigma / S_k;
            sumPtF_av += nF_k * sqrt(sqrt(eta_k));
            sumPtB_av += nF_k * sqrt(sqrt(eta_k));
            sumPtF_disp += nF_k * sqrt(eta_k);
            sumPtB_disp += nB_k * sqrt(eta_k);

            delete borderCoordAll;
        }
    }
    pF_i_av = sumPtF_av / nF_i;
    pB_i_av = sumPtB_av / nB_i;
    pF_i_disp = sumPtF_disp / (nF_i * nF_i);
    pB_i_disp = sumPtB_disp / (nB_i * nB_i);
    normal_distribution<float> disNorm_ptptF(pF_i_av, ptGammaSquared*sqrt(pF_i_disp));
    normal_distribution<float> disNorm_ptptB(pB_i_av, ptGammaSquared*sqrt(pB_i_disp));
    pF_i = disNorm_ptptF(gen);
    pB_i = disNorm_ptptB(gen);

    v_nF_i.push_back(nF_i);
    v_nB_i.push_back(nB_i);
    v_pF_i.push_back(pF_i);
    v_pB_i.push_back(pB_i);
}

void Simulation::f_bCalc()
{
    //nn
    // nFnB_Av
    for(uint i = 0; i < v_nF_i.size(); i++)
    {
        nFnB_Av += (float) (v_nF_i[i]*v_nB_i[i]);
    }
    nFnB_Av /= numSim;

    // nF_Av
    for(uint i = 0; i < v_nF_i.size(); i++)
    {
        nF_Av += (float)(v_nF_i[i]);
    }
    nF_Av /= numSim;

    // nB_Av
    for(uint i = 0; i < v_nB_i.size(); i++)
    {
        nB_Av += (float)(v_nB_i[i]);
    }
    nB_Av /= numSim;

    // nF2_Av
    for(uint i = 0; i < v_nF_i.size(); i++)
    {
        nF2_Av += (float)(v_nF_i[i]*v_nF_i[i]);
    }
    nF2_Av /= numSim;

    // nF_Av2
    nF_Av2 = (float)(nF_Av*nF_Av);

    // ptpt
    // pFpB_Av
    for(uint i = 0; i < v_pF_i.size(); i++)
    {
        pFpB_Av += (float) (v_pF_i[i]*v_pB_i[i]);
    }
    pFpB_Av /= numSim;

    // pF_Av
    for(uint i = 0; i < v_pF_i.size(); i++)
    {
        pF_Av += (float)(v_pF_i[i]);
    }
    pF_Av /= numSim;

    // pB_Av
    for(uint i = 0; i < v_pB_i.size(); i++)
    {
        pB_Av += (float)(v_pB_i[i]);
    }
    pB_Av /= numSim;

    // pF2_Av
    for(uint i = 0; i < v_pF_i.size(); i++)
    {
        pF2_Av += (float)(v_pF_i[i]*v_pF_i[i]);
    }
    pF2_Av /= numSim;

    // pF_Av2
    pF_Av2 = (float)(pF_Av*pF_Av);

    // b
    b_nn = (nFnB_Av - nF_Av*nB_Av) / (nF2_Av - nF_Av2);
    b_pp = (pFpB_Av - pF_Av*pB_Av) / (pF2_Av - pF_Av2);

    data_b << b_nn << "\t" << b_pp << "\t";
}

/// data management

void Simulation::f_WriteData()
{
    data_nF_i << nF_i << "\t";
    data_nB_i << nB_i << "\t";

    data_pF_i << pF_i << "\t";
    data_pB_i << pB_i << "\t";

    /*
    ofstream TEST;
    TEST.open("TEST.txt", ios_base::app);
    for(int i = 0; i < v_x.size(); i++)
    {
        TEST << v_x[i] << "\t";
    }
    TEST << endl;

    TEST.close();
    */
}
