
#include "../stdafx.h"

using namespace std;

void dfs(const Graph & graph , Vertex * source , std::vector<Vertex *>& path);
void dfsAux(const Graph & graph , Vertex * current , std::vector<Vertex *> &path);

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


std::vector<Vertex *> findArticulationPoints(Graph * g );
vector<vector<Vertex *>> sccTarjan(Graph* g);

void firstDFSKosarajuSharir(Vertex *v, stack<Vertex *> *vertexStack) ;
void secondDFSKosarajuSharir(Vertex *v, vector<Vertex *>& res) ;
void aux_reverseGraphEdges(Graph * g) ;
vector<vector< Vertex *>> SCCkosaraju(Graph * g) ;
