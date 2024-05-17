#include <random>
#include <iomanip>
#include "Algorithms.h"
#include "algorithm"
#include "set"

#include "../data_structures/MutablePriorityQueue.h"
#include <locale>
using namespace  std;


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

void dfs(Graph & graph , Vertex * source , std::vector<Vertex *>& path){
    std::unordered_map<int, Vertex *> vertexSet = graph.getVertexSet();
    for(std::pair<int , Vertex *> pair : vertexSet){
        Vertex * v = pair.second;
        v->setVisited(false);
    }

    dfsAux(graph , source , path);
}

std::vector<Edge * > NearestNeighbour(Graph * graph , Vertex * startVertex) {

    if(!graph->makeFullyConnected()){
        std::cout << "Not possible to make fully connected\n";
        return {};
    }

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

    Vertex * curr = startVertex;
    curr->setVisited(true);

    while (true){
        std::vector<Edge *> edgesV;
        for(Edge * e : curr->getAdj()){
            edgesV.push_back(e);
        }

        std::sort(edgesV.begin(), edgesV.end(), [](const Edge * e1, const Edge* e2)
        {
            return e1->getWeight() < e2->getWeight();
        });

        bool found = false;

        for(Edge * e : edgesV){
            if(!e->getDestination()->isVisited()){
                optimalPath.push_back(e);
                e->getDestination()->setVisited(true);
                curr = e->getDestination();
                found = true;
                break;
            }
        }

        if(!found){
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

Edge * findEdgeTo(Vertex * source , Vertex * dest){
    for(Edge * e : source->getAdj()){
        if(e->getDestination() == dest){
            return e;
        }
    }
    return nullptr;
}

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

void preOrder(const Graph * graph , Vertex * source , std::vector<Vertex *>& path){
    std::unordered_map<int, Vertex *> vertexSet = graph->getVertexSet();
    for(std::pair<int , Vertex *> pair : vertexSet){
        Vertex * v = pair.second;
        v->setVisited(false);
    }

    preOrderDFSAux(graph , source , path);
}

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
        if (curr->isVisited()) continue;
        curr->setVisited(true);
        path.push_back(curr);
        // Update distances and paths to adjacent vertices

        for(Edge * e : curr->getAdj()){
            Vertex * neighbor = e->getDestination();

            if(neighbor->isVisited()){
                continue;
            }else if(neighbor->getDist() == DBL_MAX){
                e->setSelected(true);
                neighbor->setPath(e);
                neighbor->setDist(e->getWeight());
                vertexQueue.insert(neighbor);
            } else if(e->getWeight() < neighbor->getDist()){
                e->setSelected(true);
                neighbor->getPath()->setSelected(false);
                neighbor->setPath(e);
                neighbor->setDist(e->getWeight());
                vertexQueue.decreaseKey(neighbor);
            }
        }
    }

    clock1.elapsed();

    std::cout << "Prim took " << time << "second" << std::endl;
    return path;
}



/* A heurística de "Vizinho Mais Próximo" (N.N - Nearest Neighbor) é
   um algoritmo de busca heurística utilizado para encontrar um caminho aproximadamente mais curto em um grafo, especialmente em problemas de TSP (Traveling Salesman Problem). A ideia básica do algoritmo é começar de um nó inicial e, em cada etapa, selecionar o nó mais próximo que ainda não foi visitado.

    Aqui está uma descrição geral do algoritmo:

    Escolha um nó inicial como o nó atual.
    Enquanto houver nós não visitados:
    a. Encontre o vizinho mais próximo do nó atual que ainda não foi visitado.
    b. Adicione esse vizinho ao caminho.
    c. Marque o vizinho como visitado e faça dele o nó atual.
    Retorne ao nó inicial para completar o ciclo.
 */

/*Triangulolar inequality : the least distance path to reach a vertex j from i is always to reach j directly from i
   rather than go through some vertex k.
   dist(i, j ) <= dist(i,k) + dist(k , j)
   1. Algorithm: select a root vertex
   2. Find a minimum spanning tree
   3. Do preorder walk of T. and return Hamilton Cycle
 *
 *
*/

std::vector<Edge *> TriangularApproximationHeuristic(Graph * graph , Vertex * source , Vertex * dest){

    if(!graph->makeFullyConnected()){
        std::cout << "Not possible to make fully connected\n";
        return {};
    }


    std::vector<Vertex * > path = PrimMst(graph , source);

    std::vector<Vertex *> optimalRoute;
    std::vector<Edge *> bestPath;


    if(path.size() != graph->getVertexSet().size()){
        std::cout << "Not all nodes visited. There is no Hamilton";
        return {};
    }

   preOrder( graph , source, optimalRoute );

    for(int i = 0 ; i < optimalRoute.size() -1 ; i++){
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
            std::cout << (*it)->getId() << " - ";
        }
        std::cout << (*it)->getId();
    }

    auto duration_seconds = std::chrono::duration_cast<std::chrono::seconds>(end_time - start_time);
    auto duration_ms = std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time);

    std::cout << std::endl << std::endl <<  "The execution time was: " << duration_seconds.count() << " seconds" <<  std::endl;
    std::cout << "The execution time was: " << duration_ms.count() << " milliseconds" <<  std::endl;
}

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
/*
double twoOpt(std::vector<Vertex*>& path, Graph* g, int max_iterations) {
    int n = path.size();
    bool improvement = true;
    double total_cost = 0.0;
    int iteration_count = 0;

    while (improvement && iteration_count < max_iterations) {
        improvement = false;
        for (int i = 0; i < n - 2; ++i) {
            for (int j = i + 2; j < n; ++j) {
                Edge* edge_1 = g->findEdge(path[i]->getId(), path[j]->getId());
                Edge* edge_2 = g->findEdge(path[i + 1]->getId(), path[(j + 1) % n]->getId());
                Edge* edge_3 = g->findEdge(path[i]->getId(), path[i + 1]->getId());
                Edge* edge_4 = g->findEdge(path[j]->getId(), path[(j + 1) % n]->getId());

                if (edge_1 && edge_2 && edge_3 && edge_4) {
                    double delta = (edge_1->getWeight() + edge_2->getWeight()) - (edge_3->getWeight() + edge_4->getWeight());

                    if (delta < 0) {
                        std::reverse(path.begin() + i + 1, path.begin() + j + 1);
                        improvement = true;
                        break;
                    }
                }
            }
            if (improvement) break;
        }
        iteration_count++;
    }

    for (int i = 0; i < n - 1; ++i) {
        int edge_cost = g->findEdge(path[i]->getId(), path[i + 1]->getId())->getWeight();
        std::cout << path[i]->getId() << " - " << path[i + 1]->getId() << " : " << edge_cost << std::endl;
        total_cost += edge_cost;
    }


    total_cost += g->findEdge(path[n - 1]->getId(), path[0]->getId())->getWeight();
    std::cout << path[n - 1]->getId() << " - " << path[0]->getId() << " : " << g->findEdge(path[n - 1]->getId(), path[0]->getId())->getWeight() << std::endl;

    return total_cost;
}*/


void twoOptSwap(std::vector<Vertex *> &path, int i, int k) {
    std::reverse(path.begin() + i, path.begin() + k + 1);
}

double calculateCostDifference(Vertex *a, Vertex *b, Vertex *c, Vertex *d, std::vector<std::vector<double>> &edgeMatrix) {
    return edgeMatrix[a->getId()][c->getId()] + edgeMatrix[b->getId()][d->getId()] - edgeMatrix[a->getId()][b->getId()] - edgeMatrix[c->getId()][d->getId()];
}

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
                if(k == n - 1 || i == 0) continue; // Skip the first and last vertices (they are the same in a TSP tour
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

void orderEdges(Vertex * v) {

    auto comp = [] (Edge * e1, Edge * e2){
        if (e1->getDestination()->getDegree() < e2->getDestination()->getDegree()) {
            return true;
        }
        else if (e1->getDestination()->getDegree() == e2->getDestination()->getDegree() and e1->getWeight() < e2->getWeight()){
            return true;
        }
        return false;
    };

    std::vector<Edge *> edges = v->getAdj();

    std::sort(edges.begin(), edges.end(), comp);

    v->setAdj(edges);
}


bool nn_backtracking(int & size , Vertex * s, Vertex * d, std::vector<Vertex *>& path) {
    s->setVisited(true);
    path.push_back(s);
    for (Edge * e : s->getAdj()){
        if (!e->getDestination()->isVisited()){
            if (nn_backtracking(size , e->getDestination(), d, path)) return true;
        }
        else if (path.size() == size){ /// Check if the vertex is connected to the destination

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

bool nn_with_backtracking(Graph * g , Vertex * s, std::vector<Vertex * > &hamiltonian) {

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

    hamiltonian.clear();
    std::cout << "Enter\n";
    int size = g->getVertexSet().size();
    if (!nn_backtracking( size ,  s, s, hamiltonian))
        return false;

    std::cout << "Left\n";

    std::vector<std::vector<double>> edgeMatrix(size, std::vector<double>(size, DBL_MAX));

    for(auto v : g->getVertexSet()){
        v.second->setVisited(false);
        for(auto e : v.second->getAdj()){
            edgeMatrix[v.second->getId()][e->getDestination()->getId()] = e->getWeight();
        }
    }

    //twoOpt(hamiltonian , 10 , edgeMatrix);

    std::set<uint32_t> set;
    double cost = 0;
    Vertex * prev = nullptr;
    for (Vertex * x : hamiltonian){
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


bool nn_with_backtrackingAndTwoOpt(Graph * g , Vertex * s, std::vector<Vertex * > &hamiltonian){
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

    hamiltonian.clear();
    std::cout << "Enter\n";
    int size = g->getVertexSet().size();
    if (!nn_backtracking( size ,  s, s, hamiltonian))
        return false;

    std::cout << "Left\n";

    std::vector<std::vector<double>> edgeMatrix(size, std::vector<double>(size, DBL_MAX));

    for(auto v : g->getVertexSet()){
        v.second->setVisited(false);
        for(auto e : v.second->getAdj()){
            edgeMatrix[v.second->getId()][e->getDestination()->getId()] = e->getWeight();
        }
    }

    twoOpt(hamiltonian , 10 , edgeMatrix);

    std::set<uint32_t> set;
    double cost = 0;
    Vertex * prev = nullptr;
    for (Vertex * x : hamiltonian){
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
//11982915

/*
void simulatedAnnealing(Graph* g, double initial_temperature, double cooling_rate, int max_iterations) {

    auto start_time = std::chrono::steady_clock::now();

    std::vector<Vertex*> best_path;
    double best_cost = std::numeric_limits<double>::max();

    std::vector<Vertex*> current_path;
    double current_cost = 0;

    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<> dis(0, 1);

    for (int i = 0; i < max_iterations; ++i) {

        current_path.clear();
        current_cost = 0;

        std::unordered_set<Vertex*> visited_vertices;

        Vertex* current_vertex = g->findVertex(0); // Assuming it starts from vertex 0, can be adjusted later
        visited_vertices.insert(current_vertex);

        // Create a vector of edges that are not selected and exist
        while (visited_vertices.size() < g->getNumVertex()) {
            std::vector<Edge*> unvisited_edges;

            for (Edge *e : current_vertex->getAdj()) {
                if (!e->isSelected() && visited_vertices.find(e->getDestination()) == visited_vertices.end()) {
                    unvisited_edges.push_back(e);
                }
            }

            if (unvisited_edges.empty()) {
                break;
            }

            // Uniform distribution to choose a random edge
            std::uniform_int_distribution<> distributionUnvisited(0, unvisited_edges.size() - 1);
            int random_index = distributionUnvisited(gen);
            Edge* selected_edge = unvisited_edges[random_index];
            selected_edge->setSelected(true);

            current_path.push_back(current_vertex);
            current_cost += selected_edge->getWeight();

            current_vertex = selected_edge->getDestination();
            visited_vertices.insert(current_vertex);
        }

        // Check if the last vertex has a connection to the first one
        Edge* last_to_first_edge = findEdgeTo(current_vertex, g->findVertex(0));
        if (last_to_first_edge != nullptr) {
            // Add the last edge to return to the starting vertex
            last_to_first_edge->setSelected(true);
            current_path.push_back(current_vertex);
            current_path.push_back(g->findVertex(0)); // Assuming it starts from vertex 0, can be adjusted later
            current_cost += last_to_first_edge->getWeight();
        }

        // Decide whether to accept or reject the solution
        if (current_cost < best_cost && current_path.size() == g->getVertexSet().size() + 1) {
            best_path = current_path;
            best_cost = current_cost;
        }

        // Reset visited and selected status for the next iteration
        for(auto pair: g->getEdgeSet()){
            pair.second->setSelected(false);
        }
    }

   best_cost = twoOpt(best_path, g,max_iterations);

    auto end_time = std::chrono::steady_clock::now();
    auto elapsed_time = std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time).count();

    std::cout << "Best cost: " << best_cost << std::endl;
    std::cout << "Best path: ";
    for (Vertex* vertex : best_path) {
        std::cout << vertex->getId() << " ";
    }
    std::cout << std::endl;
    std::cout << "Execution time: " << elapsed_time << " milliseconds" << std::endl;
}*/


/*
    Initialization: We start by placing the ants on the starting point.
    Movement: Each ant selects a point to move to based on a probabivectoric function that takes into account the pheromone level and the heuristic information. The heuristic information can be thought of as a measure of how good a particular point is.
    Updating pheromone trails: The pheromone trail on each edge is updated based on the quality of the solution found by the ant that traversed that edge.
    Termination: We stop the algorithm after a certain number of iterations or when a satisfactory solution is found.

std::vector<Edge *> ACO_TSP(Graph *graph, Vertex *startVertex, int numAnts, double evaporationRate, double alpha, double beta, int maxIterations, double Q, double elitistRatio) {

    std::vector<Edge *> bestTour;
    double bestTourLength = std::numeric_limits<double>::max();
    std::unordered_map<Edge *, double> probabilityTable;

    // Initialization: Set initial pheromone levels based on edge weights
    for (auto &pair : graph->getEdgeSet()) {
        pair.second->setPheromones(1.0 / pair.second->getWeight());
        probabilityTable[pair.second] = 1.0;
    }

    for (int iter = 0; iter < maxIterations; ++iter) {
        for (int antIndex = 0; antIndex < numAnts; ++antIndex) {
            Vertex *currentVertex = startVertex;
            startVertex->setVisited(true);
            std::vector<Edge *> tour;
            double tourLength = 0.0;
            int tries = graph->getNumVertex() ;
            // Move ant to the next city
            while (tour.size() < graph->getNumVertex()-1 && tries > 0) {

                // Calculate probabilities for adjacent cities
                double sumOfPheromoneLevels = 0.0;

                for (Edge * e : currentVertex->getAdj()) {
                    Vertex *nextVertex = e->getDestination();
                    if (!nextVertex->isVisited()) {
                        double pheromoneLevel = e->getPheromones();
                        double visibility = 1.0 / e->getWeight(); // Inverse of distance
                        double probability = pow(pheromoneLevel, alpha) * pow(visibility, beta);

                        probabilityTable[e] = probability;
                        sumOfPheromoneLevels += probability;
                    }
                }

                // Select next city based on probabilities
                Edge *selectedEdge = nullptr;
                double randValue = static_cast<double>(std::rand()) / RAND_MAX;
                double cumulativeProbability = 0.0;

                for (Edge *e : currentVertex->getAdj()) {
                    if (!e->getDestination()->isVisited()) {
                        cumulativeProbability += probabilityTable[e] / sumOfPheromoneLevels;
                        if (randValue <= cumulativeProbability) {
                            selectedEdge = e;
                            break;
                        }
                    }
                }

                if (selectedEdge) {
                    tour.push_back(selectedEdge);
                    tourLength += selectedEdge->getWeight();
                    selectedEdge->getDestination()->setVisited(true);
                    currentVertex = selectedEdge->getDestination();
                }

                tries--;
                //std::cout << tries << std::endl;
            }

            Edge * lastEdge = findEdgeTo(currentVertex , startVertex);

            if(lastEdge == nullptr){
                double distance = Harverstein(currentVertex->getLongitude() , currentVertex->getLatitude() , startVertex->getLongitude() , startVertex->getLatitude());
                lastEdge = new Edge(currentVertex , startVertex , -1 , distance);
            }

            tour.push_back(lastEdge);
            // Update pheromone trails along the tour path (local pheromone update)
            for (Edge *edge : tour) {
                double newPheromoneLevel = (1.0 - evaporationRate) * edge->getPheromones() + Q / tourLength;
                edge->setPheromones(newPheromoneLevel);
            }

            // Update best tour
            if (tourLength < bestTourLength) {
                bestTour = tour;
                bestTourLength = tourLength;
            }

            // Reset visited status for next iteration
            for (auto &pair : graph->getVertexSet()) {
                pair.second->setVisited(false);
            }
        }

        // Update global pheromone trails using elitist strategy
        for (Edge *edge : bestTour) {
            double newPheromoneLevel = (1.0 - evaporationRate) * edge->getPheromones() + elitistRatio * Q / bestTourLength;
            edge->setPheromones(newPheromoneLevel);
        }
    }

    return bestTour;
}
 */