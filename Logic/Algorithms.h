
#include "../stdafx.h"

std::vector<Edge * > NearestNeighbour(Graph * graph , Vertex* startVertex);

double Harverstein(double longitude1, double latitude1, double longitude2, double latitude2);

Vertex * minKey(const Graph * graph ,  const std::vector<double>& key, const std::unordered_set<Vertex *>& mstSet);
std::vector<Vertex*> PrimMst(const Graph* graph , Vertex * sourceVertex);

std::vector<Edge *> TriangularApproximationHeuristic(Graph * graph , Vertex * source , Vertex * dest);

std::vector<Edge *> ACO_TSP(Graph *graph, Vertex *startVertex, int numAnts, double evaporationRate, double alpha, double beta, int maxIterations, double Q, double elitistRatio);

Edge * findEdgeTo(Vertex * source , Vertex * dest);

std::vector<Edge*> larkeWrightSavings(Graph * graph , Vertex * sourceVertex);

std::vector<Edge *> ChristofidesAlgo(Graph * g , Vertex * source);

void tspBruteForce(Graph* g);

void tspBacktrackingBruteForce(Graph* g,Vertex* curr,double curr_cost,int n_visited,double&min_cost,std::vector<Vertex*> &path);

std::vector<Vertex*> linKernighan(Graph& graph);