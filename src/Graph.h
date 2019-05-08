#ifndef GRAPH_H
#define GRAPH_H
#include<iostream> 
#include <list> 
using namespace std; 
  

// Graph class represents a undirected graph 
// using adjacency list representation 
class Graph 
{ 
    public: 

        Graph(int V);   // Constructor 
        void addEdge(int v, int w); 
        void connectedComponents(); 

    private:

        int V; // No. of vertices 
        list<int> *adj; // Pointer to an array containing adjacency lists     
        void DFSUtil(int v, bool visited[]); // A function used by DFS 
}; 


#endif // GRAPH_H