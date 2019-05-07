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

        Simulation(std::mt19937, float, uint);
        virtual ~Simulation();
        void f_nGen();
        void f_GenerateXY(); // for uniform
        void f_FillGraph();
        void f_FindConnComp();
        void f_dfs(usint);
        void f_FindMulPtFB(uint);
        void f_sCalc();
        
        inline usint f_in_PDF_N(usint);
        inline float f_in_distXY(float, float, float, float, usint);
        inline float f_in_PDF(float, float);

    private:

        std::mt19937 gen;

        constexpr static float R = 7.5, rs = 0.225; // nucleus' radius, string's radius (def. 0.03)
        constexpr static float S_0 = M_PI*R*R, stringSigma = M_PI*rs*rs; // nucleus' and string's area
        constexpr static float ptGammaSquared = 0.5; // proportionality coefficient in pt distribution

        usint N_mean, N;
        float eta_mean;
        uint n_sim;
        std::vector<float> v_x; // vector of x coordinates
        std::vector<float> v_y; // vector of y coordinate
        std::vector<std::vector<usint>> g_connGraph; // graph for conn comps
        std::vector<bool> v_used; // vector of used comps in dfs
        std::vector<usint> v_comp; // vector for temp conn comps
        std::vector<std::vector<usint>> v_compData; // matrix for conn comps

        std::vector<float> v_S; // vector of areas
        float S_Av, S2_Av, S_Disp, S_Omega;

        std::ofstream data_S;

        clock_t tLog;
};

#endif // SIMULATION_H
