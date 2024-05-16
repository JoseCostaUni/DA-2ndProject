
#include "../stdafx.h"

using namespace std;

struct WeightedEdge {
    Edge* edge;
    double distance;

    WeightedEdge(Edge* e, double d) : edge(e), distance(d) {}

    // Comparator for priority queue
    bool operator>(const WeightedEdge& other) const {
        return distance > other.distance;
    }
};


struct CompareEdgeWeight {
public:
    bool operator()(Edge * e1, Edge * e2) {
        return e1->getWeight() > e2->getWeight();
    }
};

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

void simulatedAnnealing(Graph* g, double initial_temperature, double cooling_rate, int max_iterations);

void XNearestNeighbor(Graph * g, Vertex* src, int x);

void orderEdges(Vertex * v);
bool nn_with_backtracking(Graph * g , Vertex * s, std::vector<Vertex * > &hamiltonian) ;
bool nn_backtracking(Vertex * s, Vertex * d, std::vector<Vertex *>& path);