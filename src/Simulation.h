#ifndef SIMULATION_H
#define SIMULATION_H
#include <fstream>
#include <random>
#include <vector>

#define _USE_MATH_DEFINES
#include <cmath>
    #define M_PI       3.14159265358979323846  
//#include <malloc.h>

using uint = unsigned int;
using sint = short int;
using usint = unsigned short int;


class Simulation
{
    public:

        Simulation(std::mt19937, uint);
        virtual ~Simulation();
        void f_nGen();
        void f_GenerateXY();
        void f_FillGraph();
        void f_FindConnComp();
        void f_dfs(usint);

        
        inline usint f_in_PDF_N(usint);
        inline float f_in_distXY(float, float, float, float, usint);
        inline float f_in_PDF(float, float);

    private:

        usint N_mean, N;
        std::vector<float> v_x = std::vector<float>(N);
        std::vector<float> v_y = std::vector<float>(N);
        std::vector<std::vector<usint>> g_connGraph = std::vector<std::vector<usint>(N)>(N);
        std::vector<bool> v_used = std::vector<bool>(N); // vector of used comps in dfs
        
        std::mt19937 gen;
};

#endif // SIMULATION_H
