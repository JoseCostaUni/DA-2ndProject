
#include "../stdafx.h"

std::vector<Vertex * > NearestNeighbour(Graph * graph , Vertex* startVertex);

Vertex * minKey(const Graph * graph ,  const std::vector<double>& key, const std::unordered_set<Vertex *>& mstSet);
std::vector<Vertex*> PrimMst(const Graph* graph , Vertex * sourceVertex);

std::vector<Edge *> TriangularApproximationHeuristic(Graph * graph , Vertex * source , Vertex * dest);

void tspBruteForce(Graph* g);

void tspBacktrackingBruteForce(Graph* g,Vertex* curr,double curr_cost,int n_visited,double&min_cost,std::vector<Vertex*> &path);