#include "Simulation.h"
#include <iostream>
#include <random>
#include <algorithm>
#include <time.h>
#include <omp.h>

using namespace std;

Simulation::Simulation(int nSim)
{
    numSim = nSim + 1;
    cout << "Simulation: " << numSim << endl;

    data_Mul.open("data_Mul.txt", ios_base::app);
    data_Mul << endl << numSim << "\t\t";

    //data_Mul0.open("data_Mul0.txt", ios_base::app);
    //data_Mul0 << endl << numSim << "\t\t";

    data_MulRatio.open("data_MulRatio.txt", ios_base::app);
    data_MulRatio << endl << numSim << "\t\t";

    data_Eta.open("data_Eta.txt", ios_base::app);
    data_Eta << endl << numSim << "\t\t";

    for(int i = 0; i < a_nRange.size(); i++)
    {
        N = a_nRange[i];
        f_GenerateXY_Neym();
        f_FillGraph();
        f_FindConnComp();
        f_FindMulRatio();
        f_WriteData();
    }

    data_Mul.close();
    //data_Mul0.close();
    data_MulRatio.close();
    data_Eta.close();
}

Simulation::~Simulation() {}

/// initialization

void Simulation::f_GenerateXY_Neym()
{
    v_x.clear();
    v_y.clear();
    double xNeym, yNeym, z;
    int i = 0;
    mt19937 gen(time(0));
    uniform_real_distribution<double> dis1(-1.0, 1.0);
    uniform_real_distribution<double> dis2(0, 1.0);
    while(i < N)
    {
        xNeym = dis1(gen);
        yNeym = dis1(gen);
        //yNeym = 2 * (double)(rand()) / RAND_MAX - 1;
        //z = (double)rand() / RAND_MAX;
        z = dis2(gen);
        if(z < f_in_PDF(xNeym, yNeym))
        {
            //x[i] = xNeym;
            //y[i] = yNeym;
            v_x.push_back(xNeym);
            v_y.push_back(yNeym);
            ++i;
        }
    }
}

void Simulation::f_FillGraph()
{
    g_connGraph.clear();
    for(int i = 0; i < N; i++)
    {
        vector<int> temp;
        for(int j = 0; j < N; j++)
            temp.push_back(-1);
        g_connGraph.push_back(temp);
        temp.clear();
    }

    #pragma omp parallel for
    for(int i = 0; i < N; i++)
    {
        for(int j = 0; j < N; j++)
        {
            if(i != j)
                if(f_in_distXY(v_x[i], v_x[j], v_y[i], v_y[j], 2) <= 0)
                {
                    g_connGraph[i][j] = j; // the graph made out of -1's and numbers of all the centers of strings overlapped with j'th

                }
        }
    }
}

/// conn comp calc

void Simulation::f_FindConnComp() // used in constructor
{
    v_compData.clear(); // clear the vcm
    v_comp.clear(); // clear the temp cluster array
    v_used.clear(); // clear the trigger marks
    for(int i=0; i < N; ++i)
		v_used.push_back(false); // make all the coordinates unused
	for(int i=0; i < N; ++i)
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

/// dfs

void Simulation::f_dfs(int v)
{
    v_used[v] = true; // set v'th coordinate as being searched for near strings, overlapping with v'th
	v_comp.push_back(v); // push v'th string to a component being made
	for(int i = 0; i < N; ++i)
	{
	    if(g_connGraph[v][i] > 0) // do if not equal to -1
        {
            int to = g_connGraph[v][i]; // get a neighbor from the graph
            if(! v_used[to]) // check if a neighbor was not checked yet
                f_dfs(to); // recursively repeat dfs with all the found neighbors not checked yet (unused)
        }
	}
}

/// math

double Simulation::f_in_PDF(double x, double y)
{
    return (-x*x - y*y + 1);
}

double Simulation::f_in_distXY(double xi, double xj, double yi, double yj, int factor)
{
    return (sqrt((xi - xj) * (xi - xj) + (yi - yj) * (yi - yj)) - factor * rs);
}



/// cluster analysis

void Simulation::f_FindMulRatio()
{
    double *borderCoordAll; // for the coordinates of the border of a cluster area: xL xR yD yU
    vector<double> v_xObs, v_yObs;
    double overallArea, MCDist, MCPointX, MCPointY, S_k, S;
    unsigned int MCNThrown, MCNAll, N_k;

    Mul = 0;
    Mul0 = 0;
    Eta = 0;
    S = 0;

    for(int clusIter = 0; clusIter < v_compData.size(); clusIter++)
    {
        if(v_compData[clusIter].size() != 0)
        {
            N_k = v_compData[clusIter].size();
            Mul0 = N;

            v_xObs.clear();
            v_yObs.clear();
            for(int obsIter = 0; obsIter < v_compData[clusIter].size(); obsIter++)
            {
                v_xObs.push_back(v_x[v_compData[clusIter][obsIter]]);
                v_yObs.push_back(v_y[v_compData[clusIter][obsIter]]);
            }

            borderCoordAll = new double [4];
            borderCoordAll[0] = *std::max_element(v_xObs.begin(), v_xObs.end()) + rs;
            borderCoordAll[1] = *std::min_element(v_xObs.begin(), v_xObs.end()) - rs;
            borderCoordAll[2] = *std::max_element(v_yObs.begin(), v_yObs.end()) + rs;
            borderCoordAll[3] = *std::min_element(v_yObs.begin(), v_yObs.end()) - rs;

            MCPointX = borderCoordAll[1] + 0.003;
            MCPointY = borderCoordAll[3] + 0.003;
            MCNThrown = 0; // overall
            MCNAll = 0; // in the string
            MCDist = 0;

            while(MCPointY <= borderCoordAll[2] - 0.003)
            {
                while(MCPointX <= borderCoordAll[0] - 0.003)
                {
                    MCNThrown++;
                    for(int i = 0; i < v_xObs.size(); i++)
                    {
                        MCDist = f_in_distXY(MCPointX, v_xObs[i], MCPointY, v_yObs[i], 1);
                        if(MCDist <= 0)
                        {
                            MCNAll++;
                            break;
                        }
                    }
                    MCPointX = MCPointX + 0.01;
                }
                MCPointY = MCPointY + 0.01;
                MCPointX = borderCoordAll[1];
            }
            overallArea = (borderCoordAll[0] - borderCoordAll[1])*(borderCoordAll[2] - borderCoordAll[3]);

            S_k = overallArea * (1 + 0.087 / N_k) * MCNAll / MCNThrown;
            //cout << MCNAll << endl << MCNThrown << endl << S_k << endl;
            Mul = Mul + sqrt(N_k * S_k / stringSigma);
            S = S + S_k;

            delete borderCoordAll;
        }
    }
    //cout << Mul << endl;
    MulRatio = Mul / Mul0;
    Eta = N * stringSigma / S;
}

/// data management

void Simulation::f_WriteData()
{
    data_Mul << Mul << "\t";
    //data_Mul0 << Mul0 << "\t";
    data_MulRatio << MulRatio << "\t";
    data_Eta << Eta << "\t";
}
