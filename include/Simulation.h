#ifndef SIMULATION_H
#define SIMULATION_H
#include <vector>
#include <array>
#include <fstream>
#include <random>

#define _USE_MATH_DEFINES
#include <cmath>
#define M_PI       3.14159265358979323846
//#include <malloc.h>

using namespace std;
using uint = unsigned int;
using sint = short int;
using usint = unsigned short int;

class Simulation
{
    public:

        Simulation(mt19937, uint);  //void f_nRangeFill();
        virtual ~Simulation();
        void f_nGen();
        void f_GenerateXY(); // for uniform
        void f_GenerateXY_Neym();
        void f_FillGraph();
        void f_FindConnComp();
        void f_FindMulPtFB();
        void f_WriteData();
        void f_bCalc();

        void f_dfs(usint v);
        usint f_in_PDF_N(usint);
        float f_in_PDF(float, float);
        //double f_in_PDF_mul(double);
        float f_in_distXY(float, float, float, float, usint);
        // void f_NeymPDF();

        mt19937 gen;

    private:

        //unsigned int numString; // number of strings per one series of simulations
        float x0, y0; // x, y offset
        constexpr static float rs = 0.03, stringSigma = 3.14*rs*rs; // string's radius (def. 0.03) and area
        constexpr static float ptGammaSquared = 0.5; // proportionality coefficient in pt distribution
        //const unsigned int nRange[] = {2, 100, 500, 750, 1000, 1250, 1500, 2000, 2500, 3000, 3500, 4000, 4500, 5000, 6000, 8000, 10000, 12000}; // N range
        //unsigned int nRangeSize = sizeof(nRange)/sizeof(*nRange); // N range size

        //vector<unsigned int> v_nRange; // N range vector
        array<usint, 7> a_nRange{ {100, 250, 550, 1100, 3300, 5600, 12200} }; // N range array: rho = 3 for RHIC, 11 for LHC; N = rho / rs^2
        //array<usint, 1> a_nRange{ {100} };

        //array<usint, 1> a_nRange{ {24000} }; // FOR TESTS

        usint N, N_mean; // number of strings (randomly generated if necessary), mean for N generator
        //float N_disp; // N dispersion power
        //const unsigned int numSim = 5000; // number of simulations for one N
        uint numSim; // amount of simulations
        vector<float> v_x; // vector of x coordinate
        vector<float> v_y; // vector of y coordinate

        vector<vector<usint> > g_connGraph; // graph for conn comps
        //vector<int> temp; // temporary vector for the graph filling
        vector<bool> v_used; // vector of used comps in dfs
        vector<usint> v_comp; // vector for temp conn comps
        vector<vector<usint> > v_compData; // matrix for conn comps

        //double nA; // particle multiplicity if the strings overlap
        uint nF_i, nB_i; // F, B for 1 simulation
        vector<uint> v_nF_i, v_nB_i; // collecting nF and nB

        float sumPtF_av, sumPtB_av, sumPtF_disp, sumPtB_disp; // sum in av & disp pt calc
        float pF_i_av, pB_i_av, pF_i_disp, pB_i_disp; // av and disp pt calc
        float pF_i, pB_i; // F, B for 1 simulation
        vector<float> v_pF_i, v_pB_i; // collecting pF and pB


        float nFnB_Av, nF_Av, nB_Av, nF2_Av, nF_Av2; // averaging for the multiplicities
        float pFpB_Av, pF_Av, pB_Av, pF2_Av, pF_Av2; // averaging for the pt
        float b_nn, b_pp; //

        ofstream data_nF_i, data_nB_i, data_pF_i, data_pB_i, data_b;

};

#endif // SIMULATION_H
