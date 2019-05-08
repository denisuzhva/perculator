#ifndef SIMULATION_H
#define SIMULATION_H
#include <random>
#include <vector>
#include <fstream>
#include <time.h>

#define _USE_MATH_DEFINES
#define M_PI       3.14159265358979323846 


using uint = unsigned int;
using sint = short int;
using usint = unsigned short int;


class Simulation
{
    public:

        Simulation(std::mt19937, usint, uint);
        virtual ~Simulation();
        void f_nGen();
        void f_GenerateXY(); // for uniform
        void f_FillGraph();
        void f_FindConnComp();
        void f_dfs(usint);
        void f_FindMulPtFB(uint);
        void f_bCalc();
        
        inline usint f_in_PDF_N(usint);
        inline float f_in_distXY(float, float, float, float, usint);
        inline float f_in_PDF(float, float);

    private:

        std::mt19937 gen;

        constexpr static float R = 7.5, rs = 0.225; // nucleus' radius, string's radius
        constexpr static float S_0 = M_PI*R*R, stringSigma = M_PI*rs*rs; // area of a nucleus and a string
        constexpr static float ptGammaSquared = 0.5; // proportionality coefficient in pt distribution

        usint N_mean, N;
        uint n_sim;
        std::vector<float> v_x; // vector of x coordinates
        std::vector<float> v_y; // vector of y coordinate
        std::vector<std::vector<usint>> g_connGraph; // graph for conn comps
        std::vector<bool> v_used; // vector of used comps in dfs
        std::vector<usint> v_comp; // vector for temp conn comps
        std::vector<std::vector<usint>> v_compData; // matrix for conn comps

        float nF_i, nB_i, pF_i, pB_i; // F, B for 1 simulation
        float sumPtF_av, sumPtB_av, sumPtF_disp, sumPtB_disp; // sum in av & disp pt calc
        float pF_i_av, pB_i_av, pF_i_disp, pB_i_disp; // av and disp pt calc

    	std::vector<uint> nF_vec, nB_vec;
    	std::vector<float> pF_vec, pB_vec;
        
        float nFnB_av, nF_av, nB_av, nF2_av, nF_av2; // averaging for the multiplicities
        float pFpB_av, pF_av, pB_av, pF2_av, pF_av2; // averaging for the pt
        float b_nn, b_pp; 

        std::ofstream data_b;

        clock_t tLog;
};

#endif // SIMULATION_H
