#ifndef SIMULATION_H
#define SIMULATION_H
#include <vector>
#include <array>
#include <fstream>
using namespace std;

class Simulation
{
    public:

        Simulation(int);
        virtual ~Simulation();
        //void f_nRangeFill();
        void f_GenerateXY_Neym();
        void f_FillGraph();
        void f_FindConnComp();
        void f_FindMulRatio();
        void f_WriteData();

        void f_dfs(int v);
        double f_in_PDF(double, double);
        double f_in_distXY(double, double, double, double, int);
        // void f_NeymPDF();

    private:

        unsigned int numSim; // number of a simulation
        double x0, y0; // x, y offset
        constexpr static double rs = 0.03, stringSigma = 3.14*rs*rs; // string's radius (def. 0.03) and area
        //const unsigned int nRange[] = {2, 100, 500, 750, 1000, 1250, 1500, 2000, 2500, 3000, 3500, 4000, 4500, 5000, 6000, 8000, 10000, 12000}; // N range
        //unsigned int nRangeSize = sizeof(nRange)/sizeof(*nRange); // N range size

        //vector<unsigned int> v_nRange; // N range vector
        array<int, 14> a_nRange{ {50, 100, 500, 750, 1000, 1250, 1500, 2000, 3000, 4000, 5000, 6000, 9000, 12000} }; // N range array
        //array<int, 1> a_nRange{ {1} };
        unsigned int N; // N iterating over nRange
        //double *x, *y; // arrays of coordinates
        vector<double> v_x; // vector of x coord
        vector<double> v_y; // vector of y coord

        vector<vector<int> > g_connGraph; // graph for conn comps
        vector<bool> v_used; // vector of used comps in dfs
        vector<int> v_comp; // vector for temp conn comps
        vector<vector<int> > v_compData; // matrix for conn comps
        //int maxCompIndex; // index of the max cluster
        //int maxCompSize; // size of the max cluster

        double Mul; // particle multiplicity if the strings overlap
        double Mul0; // particle multiplicity if the strings don't overlap
        double MulRatio; // ratio between Mul and Mul0
        double Eta; // string density

        ofstream data_Mul, data_MulRatio, data_Eta;
};

#endif // SIMULATION_H
