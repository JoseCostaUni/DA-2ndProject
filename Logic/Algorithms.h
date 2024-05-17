
#include "../stdafx.h"
using namespace std;

double Harverstein(double longitude1, double latitude1, double longitude2, double latitude2);
std::vector<Vertex*> PrimMst(const Graph* graph , Vertex * sourceVertex);
Edge * findEdgeTo(Vertex * source , Vertex * dest);

void dfs(const Graph & graph , Vertex * source , std::vector<Vertex *>& path);
void dfsAux(const Graph & graph , Vertex * current , std::vector<Vertex *> &path);

std::vector<Edge * > NearestNeighbour(Graph * graph , Vertex* startVertex);

void orderEdges(Vertex * v);
bool nn_with_backtracking(Graph * g , Vertex * s, std::vector<Vertex * > &hamiltonian) ;
bool nn_with_backtrackingAndTwoOpt(Graph * g , Vertex * s, std::vector<Vertex * > &hamiltonian);
bool nn_backtracking(int & size , Vertex * s, Vertex * d, std::vector<Vertex *>& path);


void tspBruteForce(Graph* g);
void tspBacktrackingBruteForce(Graph* g,Vertex* curr,double curr_cost,int n_visited,double&min_cost,std::vector<Vertex*> &path);

std::vector<Edge *> TriangularApproximationHeuristic(Graph * graph , Vertex * source , Vertex * dest);

std::vector<Edge *> ACO_TSP(Graph *graph, Vertex *startVertex, int numAnts, double evaporationRate, double alpha, double beta, int maxIterations, double Q, double elitistRatio);
