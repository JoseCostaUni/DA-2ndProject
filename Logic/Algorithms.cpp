#include <random>
#include <iomanip>
#include "Algorithms.h"
#include "algorithm"
#include "set"

#include "../data_structures/MutablePriorityQueue.h"
#include <locale>
using namespace  std;



/**
 * @brief Auxiliary function for Depth-First Search (DFS) traversal.
 * @param graph Reference to the graph.
 * @param current Pointer to the current vertex.
 * @param path Reference to the vector storing the DFS path.
 * @complexity O(V + E), where V is the number of vertices and E is the number of edges.
 */
void dfsAux(Graph & graph , Vertex * current , std::vector<Vertex *> &path){

    path.push_back(current);
    current->setVisited(true);

    for(Edge * e : current->getAdj()){

        if(!e->getDestination()->isVisited()){
            Vertex *nextVertex = e->getDestination();
            dfsAux(graph, nextVertex, path);
        }
    }

    //current->setVisited(false); dont know if i have to put
}

/**
 * @brief Performs Depth-First Search (DFS) on the graph.
 * @param graph Reference to the graph.
 * @param source Pointer to the source vertex.
 * @param path Reference to the vector storing the DFS path.
 * @complexity O(V + E), where V is the number of vertices and E is the number of edges.
 */

void dfs(Graph & graph , Vertex * source , std::vector<Vertex *>& path){
    std::unordered_map<int, Vertex *> vertexSet = graph.getVertexSet();
    for(std::pair<int , Vertex *> pair : vertexSet){
        Vertex * v = pair.second;
        v->setVisited(false);
    }

    dfsAux(graph , source , path);
}


/**
 * @brief Implements the Nearest Neighbour heuristic for TSP.
 * @param graph Pointer to the graph.
 * @param startVertex Pointer to the starting vertex.
 * @return Vector of edges representing the approximate optimal path.
 * @complexity O(V^2 log V), where V is the number of vertices.
 */
std::vector<Edge * > NearestNeighbour(Graph * graph , Vertex * startVertex) {

    std::unordered_map<int, Vertex *> vertexSet = graph->getVertexSet();

    std::vector<Edge * > optimalPath;

    for(std::pair<int , Vertex *> pair : vertexSet){
        Vertex * v = pair.second;
        v->setDist(DBL_MAX);
        v->setVisited(false);
        v->setPath(nullptr);
        for (Edge * e : v->getAdj()){
            e->setSelected(false);
        }
    }

    Vertex * currentVertex = startVertex;
    currentVertex->setVisited(true);

    while (true){
        std::vector<Edge *> edgesV;
        for(Edge * e : currentVertex->getAdj()){
            edgesV.push_back(e);
        }

        std::sort(edgesV.begin(), edgesV.end(), [](const Edge * e1, const Edge* e2)
        {
            return e1->getWeight() < e2->getWeight();
        });

        bool foundEdge = false;

        for(Edge * edge : edgesV){
            if(!edge->getDestination()->isVisited()){
                edge->getDestination()->setVisited(true);
                optimalPath.push_back(edge);
                currentVertex = edge->getDestination();
                foundEdge = true;
                break;
            }
        }

        if(!foundEdge){
            break;
        }
    }

    if(!optimalPath.empty()){
        Edge * last = findEdgeTo(optimalPath.back()->getDestination() , startVertex);
        optimalPath.push_back(last);
    }

    if(optimalPath.size() != graph->getVertexSet().size()){
        std::cout << "No hamiltion path";
        return {};
    }

    for(auto e : optimalPath){
        if(e == nullptr){
            std::cout << "No hamiltion path";
            return {};
        }
    }

    return optimalPath;
}

/**
 * @brief Calculates the Haversine distance between two geographic coordinates.
 * @param longitude1 Longitude of the first point.
 * @param latitude1 Latitude of the first point.
 * @param longitude2 Longitude of the second point.
 * @param latitude2 Latitude of the second point.
 * @return Distance in meters.
 * @complexity O(lon(n)) where N is the number whose square root is being calculated.
 */
double Harverstein(double longitude1, double latitude1, double longitude2, double latitude2) {
    double lon1_rad = longitude1 * M_PI / 180.0;
    double lat1_rad = latitude1 * M_PI / 180.0;
    double lon2_rad = longitude2 * M_PI / 180.0;
    double lat2_rad = latitude2 * M_PI / 180.0;

    double dlon = lon2_rad - lon1_rad;
    double dlat = lat2_rad - lat1_rad;
    double a = pow(sin(dlat / 2), 2) + cos(lat1_rad) * cos(lat2_rad) * pow(sin(dlon / 2), 2);
    double c = 2 * atan2(sqrt(a), sqrt(1 - a));

    const double earth_radius_km = 6371000.0;

    double distance = earth_radius_km * c;
    return distance;
}

/**
 * @brief Finds an edge between two vertices.
 * @param source Pointer to the source vertex.
 * @param dest Pointer to the destination vertex.
 * @return Pointer to the edge if found, otherwise nullptr.
 * @complexity O(E), where E is the number of edges of the source vertex.
 */
Edge * findEdgeTo(Vertex * source , Vertex * dest){
    for(Edge * e : source->getAdj()){
        if(e->getDestination() == dest){
            return e;
        }
    }
    return nullptr;
}

/**
 * @brief Auxiliary function for pre-order DFS traversal.
 * @param graph Pointer to the graph.
 * @param current Pointer to the current vertex.
 * @param path Reference to the vector storing the pre-order DFS path.
 * @complexity O(V + E), where V is the number of vertices and E is the number of edges.
 */
void preOrderDFSAux(const Graph * graph , Vertex * current , std::vector<Vertex *> &path){

    path.push_back(current);
    current->setVisited(true);

    for(Edge * e : current->getAdj()){
        if(!e->getDestination()->isVisited() && e->isSelected()){
            Vertex *nextVertex = e->getDestination();
            preOrderDFSAux(graph, nextVertex, path);
        }
    }

    //current->setVisited(false); dont know if i have to put
}

/**
 * @brief Performs pre-order traversal on the graph.
 * @param graph Pointer to the graph.
 * @param source Pointer to the source vertex.
 * @param path Reference to the vector storing the pre-order path.
 * @complexity O(V + E), where V is the number of vertices and E is the number of edges.
 */
void preOrder(const Graph * graph , Vertex * source , std::vector<Vertex *>& path){

    std::unordered_map<int, Vertex *> vertexSet = graph->getVertexSet();
    for(std::pair<int , Vertex *> pair : vertexSet){
        Vertex * v = pair.second;
        v->setVisited(false);
    }

    preOrderDFSAux(graph , source , path);
}


/**
 * @brief Implements Prim's algorithm to find the Minimum Spanning Tree (MST) of the graph.
 * @param graph Pointer to the graph.
 * @param sourceVertex Pointer to the source vertex.
 * @return Vector of vertices representing the MST path.
 * @complexity O(E log V), where V is the number of vertices and E is the number of edges.
 */
std::vector<Vertex * > PrimMst(const Graph* graph , Vertex * sourceVertex){

    Clock clock1;
    clock1.start();

    std::vector<Vertex *> path;

    std::unordered_map<int, Vertex *> vertexSet = graph->getVertexSet();

    for(std::pair<int , Vertex *> pair : vertexSet){
        Vertex * v = pair.second;
        v->setDist(DBL_MAX);
        v->setVisited(false);
        v->setPath(nullptr);
        for (Edge * e : v->getAdj()){
            e->setSelected(false);
        }
    }

    MutablePriorityQueue<Vertex> vertexQueue;

    sourceVertex->setDist(0);
    vertexQueue.insert(sourceVertex);

    while (!vertexQueue.empty()){
        Vertex * curr = vertexQueue.extractMin();

        if (curr->isVisited())
            continue;

        curr->setVisited(true);
        path.push_back(curr);

        for(Edge * e : curr->getAdj()){
            Vertex * neighbour = e->getDestination();

            if(neighbour->isVisited()){
                continue;
            }else if(neighbour->getDist() == DBL_MAX){
                e->setSelected(true);
                neighbour->setPath(e);
                neighbour->setDist(e->getWeight());
                vertexQueue.insert(neighbour);
            } else if(e->getWeight() < neighbour->getDist()){
                e->setSelected(true);
                neighbour->getPath()->setSelected(false);
                neighbour->setPath(e);
                neighbour->setDist(e->getWeight());
                vertexQueue.decreaseKey(neighbour);
            }
        }
    }

    clock1.elapsed();

    return path;
}

/**
 * @brief Solves the TSP using the Triangular Approximation Heuristic.
 * @param graph Pointer to the graph representing the cities and paths.
 * @param source Pointer to the source vertex (starting city).
 * @param dest Pointer to the destination vertex.
 * @return Vector of edges representing the approximate TSP path.
 * @complexity O(V^2 log V), where V is the number of vertices.
 */
std::vector<Edge *> TriangularApproximationHeuristic(Graph * graph , Vertex * source , Vertex * dest){

    std::vector<Vertex * > path = PrimMst(graph , source);

    std::vector<Vertex *> optimalRoute;
    std::vector<Edge *> bestPath;


    if(path.size() != graph->getVertexSet().size()){
        std::cout << "Not all nodes visited. There is no Hamilton";
        return {};
    }

   preOrder(graph , source, optimalRoute );

    for(int i = 0 ; i < optimalRoute.size() - 1 ; i++){

        Vertex * v = optimalRoute[i];
        Vertex * nextV = optimalRoute[i+1];

        Edge * e = findEdgeTo(v , nextV);

        if(e == nullptr){
            if(v->getLatitude() == DBL_MAX || v->getLongitude() == DBL_MAX)
                return {};

            double distance = Harverstein(v->getLongitude() , v->getLatitude() , nextV->getLongitude() , nextV->getLatitude());
            e = new Edge(v , nextV , -1 , distance);
        }
        bestPath.push_back(e);


    }
    Vertex * lastV = optimalRoute.back();
    Edge * returnEdge = findEdgeTo(lastV , source);

    if(returnEdge == nullptr){
        if(lastV->getLatitude() == DBL_MAX || lastV->getLongitude() == DBL_MAX)
            return {};

        double distance = Harverstein(lastV->getLongitude() , lastV->getLatitude() , source->getLongitude() , source->getLatitude());
        returnEdge = new Edge(lastV , source , -1 , distance);
    }

    bestPath.push_back(returnEdge);

    return bestPath;
}
/**
 * @brief Solves the TSP using brute force with backtracking.
 * @param g Pointer to the graph representing the cities and paths.
 * @complexity O(V!), where V is the number of vertices.
 */
void tspBruteForce(Graph* g){
    auto start_time = std::chrono::steady_clock::now();

    std::vector<Vertex *> path;

    for(std::pair<int , Vertex *> pair: g->getVertexSet()){
        Vertex * v = pair.second;
        v->setVisited(false);
        v->setPath(nullptr);
    }

    Vertex* start = g->findVertex(0);

    start->setVisited(true);

    double cost = std::numeric_limits<double>::max();

    tspBacktrackingBruteForce(g,start,0,1,cost,path);

    std::cout << "The minimum cost obtained with backtracking was: " << cost << std::endl;

    auto end_time = std::chrono::steady_clock::now();

    std::cout << "The path chosen was: " << std::endl;

    if (!path.empty()) {
        auto it = path.begin();
        for (; it != path.end() - 1; ++it) {
            std::cout << (*it)->getId() << "->";
        }
        std::cout << (*it)->getId();
    }

    auto duration_seconds = std::chrono::duration_cast<std::chrono::seconds>(end_time - start_time);
    auto duration_ms = std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time);

    std::cout << std::endl << std::endl <<  "The execution time was: " << duration_seconds.count() << " seconds" <<  std::endl;
    std::cout << "The execution time was: " << duration_ms.count() << " milliseconds" <<  std::endl;
}

/**
 * @brief Helper function for tspBruteForce to perform backtracking.
 * @param g Pointer to the graph representing the cities and paths.
 * @param curr Pointer to the current vertex.
 * @param curr_cost Current path cost.
 * @param n_visited Number of visited vertices.
 * @param min_cost Reference to the minimum cost found.
 * @param path Reference to the vector storing the best path.
 * @complexity O(V!), where V is the number of vertices.
 */
void tspBacktrackingBruteForce(Graph* g,Vertex* curr,double curr_cost,int n_visited,double &min_cost,std::vector<Vertex*> &path){

    bool vi = false;
    double costt = 0;

    if (n_visited == g->getNumVertex()) {
        for (Edge * e  : curr->getAdj()) {
            Vertex* ver = e->getDestination();
            if(ver->getId() == g->findVertex(0)->getId()){
                costt = e->getWeight();
                vi = true;
            }
        }

        if(!vi) return;

        double cost = curr_cost + costt;

        if (cost < min_cost) {
            min_cost = cost;
            path.clear();
            path.push_back(curr);
            for (Edge* edge = curr->getPath(); edge->getSource()->getId() != g->findVertex(0)->getId(); edge = edge->getSource()->getPath()) {
                path.push_back(edge->getSource());
            }

            path.push_back(g->findVertex(0));
            std::reverse(path.begin(),path.end());
            path.push_back(g->findVertex(0));
        }
        return;
    }

    for(Edge *e : curr->getAdj()){
        Vertex* ver = e->getDestination();
        if(!ver->isVisited()){
            ver->setVisited(true);
            ver->setPath(e);
            tspBacktrackingBruteForce(g,ver,curr_cost+e->getWeight(),n_visited+1,min_cost,path);
            ver->setVisited(false);
            ver->setPath(nullptr);
        }
    }
}

/**
 * @brief Performs a 2-opt swap on the path.
 * @param path Reference to the vector of vertices representing the path.
 * @param i Index of the first vertex.
 * @param k Index of the second vertex.
 * @complexity O(k - i), where i and k are indices in the path.
 */
void twoOptSwap(std::vector<Vertex *> &path, int i, int k) {
    std::reverse(path.begin() + i, path.begin() + k + 1);
}

/**
* @brief Calculates the cost difference for a 2-opt swap.
* @param a Pointer to the first vertex.
* @param b Pointer to the second vertex.
* @param c Pointer to the third vertex.
* @param d Pointer to the fourth vertex.
* @param edgeMatrix Reference to the 2D matrix of edge weights.
* @return The cost difference resulting from the swap.
* @complexity O(1)
*/
double calculateCostDifference(Vertex *a, Vertex *b, Vertex *c, Vertex *d, std::vector<std::vector<double>> &edgeMatrix) {
    return edgeMatrix[a->getId()][c->getId()] + edgeMatrix[b->getId()][d->getId()] - edgeMatrix[a->getId()][b->getId()] - edgeMatrix[c->getId()][d->getId()];
}

/**
 * @brief Improves the path using the 2-opt algorithm.
 * @param path Reference to the vector of vertices representing the path.
 * @param maxIterations Maximum number of iterations.
 * @param edgeMatrix Reference to the 2D matrix of edge weights.
 * @complexity O(V^2 * maxIterations), where V is the number of vertices.
 */
void twoOpt(std::vector<Vertex *> &path, int maxIterations , std::vector<std::vector<double>> &edgeMatrix) {
    int n = path.size();
    if (n < 4) return; // No improvement possible if there are fewer than 4 vertices

    bool improved = true;
    int iteration = 0;

    while (improved && iteration < maxIterations) {
        improved = false;
        ++iteration;
        double current_improvement = 0.0;
        for (int i = 0; i < n - 1; ++i) {
            for (int k = i + 1; k < n; ++k) {
                if(k == n - 1 || i == 0)
                    continue;
                double delta = calculateCostDifference(path[i], path[i + 1], path[k], path[(k + 1) % n], edgeMatrix);
                if (delta < 0) {
                    twoOptSwap(path, i + 1, k);
                    improved = true;
                    current_improvement -= delta;

                    if(current_improvement > 10000000)
                        return;
                }
            }
        }
    }
}


/**
 * @brief Orders the edges of a vertex based on destination degree and edge weight.
 * @param v Pointer to the vertex whose edges need to be ordered.
 * @complexity O(E log E), where E is the number of edges of the vertex.
 */
void orderEdges(Vertex * v) {

    std::vector<Edge *> edges = v->getAdj();

    std::sort(edges.begin(), edges.end(), [] (Edge * e1, Edge * e2){
        if (e1->getDestination()->getDegree() < e2->getDestination()->getDegree()) {
            return true;
        }
        else if (e1->getDestination()->getDegree() == e2->getDestination()->getDegree() and e1->getWeight() < e2->getWeight()){
            return true;
        }
        return false;
    });

    v->setAdj(edges);
}

/**
 * @brief Nearest neighbor with backtracking to find a Hamiltonian path.
 * @param size Reference to the size of the graph.
 * @param s Pointer to the start vertex.
 * @param d Pointer to the destination vertex.
 * @param path Reference to the vector storing the Hamiltonian path.
 * @return True if a Hamiltonian path is found, false otherwise.
 * @complexity O(V!), where V is the number of vertices.
 */
bool nn_backtracking(int & size , Vertex * s, Vertex * d, std::vector<Vertex *>& path) {

    if (path.size() == size){
        if (s->getId() == d->getId()){
            path.push_back(d);
            return true;
        }
        return false;
    }

    s->setVisited(true);
    path.push_back(s);
    for (Edge * e : s->getAdj()){
        if (!e->getDestination()->isVisited()){
            if (nn_backtracking(size , e->getDestination(), d, path)) return true;
        }
        else if (path.size() == size){
            if (e->getDestination()->getId() == d->getId()){
                path.push_back(d);
                return true;
            }
        }
    }

    s->setVisited(false);
    path.pop_back();
    return false;
}

/**
 * @brief Finds a Hamiltonian path using nearest neighbor with backtracking.
 * @param g Pointer to the graph representing the cities and paths.
 * @param s Pointer to the start vertex.
 * @param optimalRoute Reference to the vector storing the Hamiltonian path.
 * @return True if a Hamiltonian path is found, false otherwise.
 * @complexity O(V!), where V is the number of vertices.
 */
bool nn_with_backtracking(Graph * g , Vertex * s, std::vector<Vertex * > &optimalRoute) {

    std::vector<Vertex *>  mst = PrimMst(g , s);

    g->populate_in_and_out_degree();

    for (auto & pair : g->getVertexSet()){
        Vertex * v = pair.second;
        orderEdges(v);
        v->setVisited(false);
    }

    std::cout << "VertexSize : " << g->getVertexSet().size() << std::endl;

    Clock clock1 = Clock();
    clock1.start();

    optimalRoute.clear();

    int size = g->getVertexSet().size();
    if (!nn_backtracking( size ,  s, s, optimalRoute))
        return false;

    std::vector<std::vector<double>> edgeMatrix(size, std::vector<double>(size, DBL_MAX));

    for(auto v : g->getVertexSet()){
        v.second->setVisited(false);
        for(auto e : v.second->getAdj()){
            edgeMatrix[v.second->getId()][e->getDestination()->getId()] = e->getWeight();
        }
    }

    std::set<int> set;
    double cost = 0;
    Vertex * prev = nullptr;

    for (Vertex * x : optimalRoute){
        std::cout << x->getId() << " - ";
        set.insert(x->getId());
        if(prev != nullptr){
            Edge * e = findEdgeTo(prev, x);
            if(e == nullptr){
                std::cout << "null edge\n";
                return false;
            }
            cost += e->getWeight();
        }
        prev = x;
    }
    std::cout << "\nSet size: " << set.size() << "\n";
    clock1.elapsed();
    std::cout.imbue(std::locale(""));
    std::cout << "Cost: " << std::fixed << std::setprecision(2) << cost << "\n";
    return true;
}

/**
 * @brief Finds a Hamiltonian path using nearest neighbor with backtracking and 2-opt optimization.
 * @param g Pointer to the graph representing the cities and paths.
 * @param s Pointer to the start vertex.
 * @param optimalRoute Reference to the vector storing the Hamiltonian path.
 * @return True if a Hamiltonian path is found, false otherwise.
 * @complexity O(V! + V^2 * maxIterations), where V is the number of vertices.
 */
bool nn_with_backtrackingAndTwoOpt(Graph * g , Vertex * s, std::vector<Vertex * > &optimalRoute){
    std::vector<Vertex *> mst = PrimMst(g , s);

    g->populate_in_and_out_degree();

    for (auto & pair : g->getVertexSet()){
        Vertex * v = pair.second;
        orderEdges(v);
        v->setVisited(false);
    }

    Clock clock1 = Clock();
    clock1.start();

    optimalRoute.clear();

    int size = g->getVertexSet().size();
    if (!nn_backtracking(size , s , s , optimalRoute))
        return false;

    std::vector<std::vector<double>> edgeMatrix(size, std::vector<double>(size, DBL_MAX));

    for(auto v : g->getVertexSet()){
        v.second->setVisited(false);
        for(auto e : v.second->getAdj()){
            edgeMatrix[v.second->getId()][e->getDestination()->getId()] = e->getWeight();
        }
    }

    twoOpt(optimalRoute , 10 , edgeMatrix);

    std::set<int> set;

    double cost = 0;
    Vertex * prev = nullptr;

    for (Vertex * x : optimalRoute){
        std::cout << x->getId() << " - ";
        set.insert(x->getId());
        if(prev != nullptr){
            Edge * e = findEdgeTo(prev, x);
            if(e == nullptr){
                std::cout << "null edge\n";
                return false;
            }
            cost += e->getWeight();
        }
        prev = x;
    }

    clock1.elapsed();
    std::cout.imbue(std::locale(""));
    std::cout << "Cost: " << std::fixed << std::setprecision(2) << cost << "\n";
    return true;
}
