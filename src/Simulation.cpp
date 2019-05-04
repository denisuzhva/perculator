#define _USE_MATH_DEFINES

#include "Simulation.h"
#include <iostream>


Simulation::Simulation(std::mt19937 genMain, uint n_mean_main) : gen(genMain), N_mean(n_mean_main) {}


Simulation::~Simulation() {}


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
    std::cout << "Current string number:\t" << N << "\tMean:\t" << N_mean;
}


void Simulation::f_GenerateXY()
{
    float ang, rad, xUni, yUni;
    std::uniform_real_distribution<float> disRad(0.0, 1.0);
    std::uniform_real_distribution<float> disAng(0.0, 2*M_PI);

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


void Simulation::f_FillGraph()
{
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