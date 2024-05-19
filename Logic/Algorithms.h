/**
 * @file Algorithms.h
 * @brief This file contains the declarations of various graph algorithms including Depth-First Search (DFS), Nearest Neighbour heuristic for TSP, Prim's algorithm for MST, and others.
 */

#include "../stdafx.h"
using namespace std;

/**
 * @brief Calculates the Haversine distance between two geographic coordinates.
 *
 * @param longitude1 Longitude of the first point.
 * @param latitude1 Latitude of the first point.
 * @param longitude2 Longitude of the second point.
 * @param latitude2 Latitude of the second point.
 * @return Distance in meters.
 * @complexity O(lon(n)) where N is the number whose square root is being calculated.
 */
double Harverstein(double longitude1, double latitude1, double longitude2, double latitude2);

/**
 * @brief Implements Prim's algorithm to find the Minimum Spanning Tree (MST) of the graph.
 *
 * @param graph Pointer to the graph.
 * @param sourceVertex Pointer to the source vertex.
 * @return Vector of vertices representing the MST path.
 * @complexity O(E log V), where V is the number of vertices and E is the number of edges.
 */
std::vector<Vertex*> PrimMst(const Graph* graph , Vertex * sourceVertex);

/**
 * @brief Finds an edge between two vertices.
 *
 * @param source Pointer to the source vertex.
 * @param dest Pointer to the destination vertex.
 * @return Pointer to the edge if found, otherwise nullptr.
 * @complexity O(E), where E is the number of edges of the source vertex.
 */
Edge * findEdgeTo(Vertex * source , Vertex * dest);

/**
 * @brief Performs Depth-First Search (DFS) on the graph.
 *
 * @param graph Reference to the graph.
 * @param source Pointer to the source vertex.
 * @param path Reference to the vector storing the DFS path.
 * @complexity O(V + E), where V is the number of vertices and E is the number of edges.
 */
void dfs(const Graph & graph , Vertex * source , std::vector<Vertex *>& path);

/**
 * @brief Auxiliary function for Depth-First Search (DFS) traversal.
 *
 * @param graph Reference to the graph.
 * @param current Pointer to the current vertex.
 * @param path Reference to the vector storing the DFS path.
 * @complexity O(V + E), where V is the number of vertices and E is the number of edges.
 */
void dfsAux(const Graph & graph , Vertex * current , std::vector<Vertex *> &path);

/**
 * @brief Implements the Nearest Neighbour heuristic for TSP.
 *
 * @param graph Pointer to the graph.
 * @param startVertex Pointer to the starting vertex.
 * @return Vector of edges representing the approximate optimal path.
 * @complexity O(V^2 log V), where V is the number of vertices.
 */
std::vector<Edge * > NearestNeighbour(Graph * graph , Vertex* startVertex);

/**
 * @brief Performs a 2-opt swap on the path.
 *
 * @param path Reference to the vector of vertices representing the path.
 * @param i Index of the first vertex.
 * @param k Index of the second vertex.
 * @complexity O(k - i), where i and k are indices in the path.
 */
void twoOptSwap(std::vector<Vertex *> &path, int i, int k);

/**
 * @brief Calculates the cost difference for a 2-opt swap.
 *
 * @param a Pointer to the first vertex.
 * @param b Pointer to the second vertex.
 * @param c Pointer to the third vertex.
 * @param d Pointer to the fourth vertex.
 * @param edgeMatrix Reference to the 2D matrix of edge weights.
 * @return The cost difference resulting from the swap.
 * @complexity O(1)
 */
double calculateCostDifference(Vertex *a, Vertex *b, Vertex *c, Vertex *d, std::vector<std::vector<double>> &edgeMatrix);

/**
 * @brief Improves the path using the 2-opt algorithm.
 *
 * @param path Reference to the vector of vertices representing the path.
 * @param maxIterations Maximum number of iterations.
 * @param edgeMatrix Reference to the 2D matrix of edge weights.
 * @complexity O(V^2 * maxIterations), where V is the number of vertices.
 */
void twoOpt(std::vector<Vertex *> &path, int maxIterations, std::vector<std::vector<double>> &edgeMatrix);

/**
 * @brief Orders the edges of a vertex based on destination degree and edge weight.
 *
 * @param v Pointer to the vertex whose edges need to be ordered.
 * @complexity O(E log E), where E is the number of edges of the vertex.
 */
void orderEdges(Vertex * v);

/**
 * @brief Nearest neighbor with backtracking to find a Hamiltonian path.
 *
 * @param size Reference to the size of the graph.
 * @param s Pointer to the start vertex.
 * @param d Pointer to the destination vertex.
 * @param path Reference to the vector storing the Hamiltonian path.
 * @return True if a Hamiltonian path is found, false otherwise.
 * @complexity O(V!), where V is the number of vertices.
 */
bool nn_with_backtracking(Graph * g , Vertex * s, std::vector<Vertex * > &hamiltonian) ;

/**
 * @brief Finds a Hamiltonian path using nearest neighbor with backtracking and 2-opt optimization.
 * @param g Pointer to the graph representing the cities and paths.
 * @param s Pointer to the start vertex.
 * @param hamiltonian Reference to the vector storing the Hamiltonian path.
 * @return True if a Hamiltonian path is found, false otherwise.
 * @complexity O(V! + V^2 * maxIterations), where V is the number of vertices.
 */
bool nn_with_backtrackingAndTwoOpt(Graph * g , Vertex * s, std::vector<Vertex * > &hamiltonian);

/**
 * @brief Finds a Hamiltonian path using nearest neighbor with backtracking.
 *
 * @param g Pointer to the graph representing the cities and paths.
 * @param s Pointer to the start vertex.
 * @param hamiltonian Reference to the vector storing the Hamiltonian path.
 * @return True if a Hamiltonian path is found, false otherwise.
 * @complexity O(V!), where V is the number of vertices.
 */
bool nn_backtracking(int & size , Vertex * s, Vertex * d, std::vector<Vertex *>& path);



/**
 * @brief Solves the TSP using brute force with backtracking.
 *
 * @param g Pointer to the graph representing the cities and paths.
 * @complexity O(V!), where V is the number of vertices.
 */
void tspBruteForce(Graph* g);

/**
 * @brief Helper function for tspBruteForce to perform backtracking.
 *
 * @param g Pointer to the graph representing the cities and paths.
 * @param curr Pointer to the current vertex.
 * @param curr_cost Current path cost.
 * @param n_visited Number of visited vertices.
 * @param min_cost Reference to the minimum cost found.
 * @param path Reference to the vector storing the best path.
 * @complexity O(V!), where V is the number of vertices.
 */
void tspBacktrackingBruteForce(Graph* g,Vertex* curr,double curr_cost,int n_visited,double&min_cost,std::vector<Vertex*> &path);


/**
 * @brief Auxiliary function for pre-order DFS traversal.
 *
 * @param graph Pointer to the graph.
 * @param current Pointer to the current vertex.
 * @param path Reference to the vector storing the pre-order DFS path.
 * @complexity O(V + E), where V is the number of vertices and E is the number of edges.
 */
void preOrderDFSAux(const Graph *graph, Vertex *current, std::vector<Vertex *> &path);

/**
 * @brief Performs pre-order traversal on the graph.
 *
 * @param graph Pointer to the graph.
 * @param source Pointer to the source vertex.
 * @param path Reference to the vector storing the pre-order path.
 * @complexity O(V + E), where V is the number of vertices and E is the number of edges.
 */
void preOrder(const Graph *graph, Vertex *source, std::vector<Vertex *> &path);


/**
 * @brief Solves the TSP using the Triangular Approximation Heuristic.
 *
 * @param graph Pointer to the graph representing the cities and paths.
 * @param source Pointer to the source vertex (starting city).
 * @param dest Pointer to the destination vertex.
 * @return Vector of edges representing the approximate TSP path.
 * @complexity O(V^2 log V), where V is the number of vertices.
 */
std::vector<Edge *> TriangularApproximationHeuristic(Graph * graph , Vertex * source , Vertex * dest);

