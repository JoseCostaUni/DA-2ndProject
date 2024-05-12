#include <random>
#include "Algorithms.h"
#include "algorithm"

#include "../data_structures/MutablePriorityQueue.h"


std::vector<Edge * > NearestNeighbour(Graph * graph , Vertex * startVertex) {
    std::unordered_map<int, Vertex *> vertexSet = graph->getVertexSet();

    std::vector<Edge * > optimalPath;

    for(std::pair<int , Vertex *> pair : vertexSet){
        Vertex * v = pair.second;
        v->setDist(DBL_MAX);
        v->setVisited(false);
        v->setPath(nullptr);
        for (auto& pair_ : v->getAdj()){
            pair_.second->setSelected(false);
        }
    }

    Vertex * curr = startVertex;
    curr->setVisited(true);

    while (true){
        std::vector<Edge *> edgesV;
        for(auto& pair : curr->getAdj()){
            edgesV.push_back(pair.second);
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
    for(auto pair : source->getAdj()){
        Edge * e = pair.second;

        if(e->getDestination() == dest){
            return e;
        }
    }
    return nullptr;
}

void preOrderDFSAux(const Graph * graph , Vertex * current , std::vector<Vertex *> &path){

    path.push_back(current);
    current->setVisited(true);

    for(auto pair : current->getAdj()){
        Edge * e = pair.second;

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
        for (auto& pair_ : v->getAdj()){
            pair_.second->setSelected(false);
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
        std::unordered_map<int , Edge *> edgeSet = curr->getAdj();

        for(std::pair<int , Edge *> edgePair : edgeSet){
            Edge * e = edgePair.second;

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

    double time = clock1.elapsed();

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

    //graph->makeFullyConnected();

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

/*
    Initialization: We start by placing the ants on the starting point.
    Movement: Each ant selects a point to move to based on a probabilistic function that takes into account the pheromone level and the heuristic information. The heuristic information can be thought of as a measure of how good a particular point is.
    Updating pheromone trails: The pheromone trail on each edge is updated based on the quality of the solution found by the ant that traversed that edge.
    Termination: We stop the algorithm after a certain number of iterations or when a satisfactory solution is found.
*/
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

            // Move ant to the next city
            while (tour.size() < graph->getNumVertex()-1) {
                std::unordered_map<int, Edge *> adjacentEdges = currentVertex->getAdj();

                // Calculate probabilities for adjacent cities
                double sumOfPheromoneLevels = 0.0;

                for (auto &edgePair : adjacentEdges) {
                    Edge *e = edgePair.second;
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

                for (auto &edgePair : adjacentEdges) {
                    Edge *e = edgePair.second;
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


void computeMWPM(Graph *g) {
    std::vector<Edge*> matching;

    for(std::pair<int,Vertex*> pair: g->getVertexSet()){
        Vertex* v = pair.second;
        if((v->getIndegree() + v->getOutdegree()) % 2 != 0){
            for(std::pair<int,Edge*> pair_edge: v->getAdj()){
                matching.push_back(pair_edge.second);
            }
        }
    }

    std::sort(matching.begin(),matching.end(),[](const Edge* edge1,const Edge* edge2) {return edge1->getWeight() < edge2->getWeight();});

    for(Edge* e : matching){
        if(!e->getDestination()->isVisited() and !e->getSource()->isVisited()){
            e->getDestination()->setVisited(true);
            e->getSource()->setVisited(true);
            e->setSelected(true);
        }
    }

}

void findEulerianCircuit(Graph* g, std::vector<Vertex*>& visited_vertices) {
    std::stack<Vertex*> vertex_stack;
    Vertex* start = g->findVertex(0);
    vertex_stack.push(start);

    while(!vertex_stack.empty()){
        Vertex* top_vertex = vertex_stack.top();
        bool unexploredEdges = false;
        for(std::pair<int,Edge*> pair : top_vertex->getAdj() ){
            Edge* e = pair.second;
            if(e->isSelected()){
               e->setSelected(false);
               vertex_stack.push(e->getDestination());
               unexploredEdges = true;
               break;
            }
        }

        if(!unexploredEdges){
            visited_vertices.push_back(top_vertex);
            vertex_stack.pop();
        }
    }

    std::reverse(visited_vertices.begin(),visited_vertices.end());

}

std::vector<Edge *> ChristofidesAlgo(Graph * g , Vertex * source){
    std::vector<Vertex *> mstPath = PrimMst(g , source);
    std::vector<Edge *> result;

    if(mstPath.size() != g->getVertexSet().size()){
        return result;
    }

    g->populate_in_and_out_degree();
    computeMWPM(g);

    std::vector<Vertex*> path;

    findEulerianCircuit(g,path);

    for(auto v : g->getVertexSet()){
        v.second->setVisited(false);
    }

    std::unordered_set<Vertex*> vertex_dup_set;
    std::vector<Vertex*> result_path;

    for(Vertex* v : path){
        if(vertex_dup_set.find(v) == vertex_dup_set.end()){
            result_path.push_back(v);
            vertex_dup_set.insert(v);
        }
    }

    auto a = result_path.size();

    for(uint64_t i = 0; i < result_path.size() - 1 ;i++){
        Edge* e = findEdgeTo(result_path[i],result_path[i+1]);
        result.push_back(e);
    }

    result.push_back(findEdgeTo(result_path.back(),source));

    auto b = result.size();
    return result;
}




























































std::vector<Edge*> larkeWrightSavings(Graph * graph , Vertex * sourceVertex) {
    std::vector<Edge*> tourEdges;

    // Step 1: Choose a random city as the hub
    int hubId = sourceVertex->getId(); // Assuming the first vertex as hub
    std::unordered_map<int, Vertex *> vset = graph->getVertexSet();
    // Step 2: Compute savings for all non-hub city pairs
    std::vector<std::pair<double, std::pair<Vertex*, Vertex*>>> savings;
    for (auto it1 = vset.begin(); it1 != vset.end(); ++it1) {
        for (auto it2 = std::next(it1); it2 != vset.end(); ++it2) {
            Vertex* v1 = it1->second;
            Vertex* v2 = it2->second;

            double saving = Harverstein(v1->getLongitude(), v1->getLatitude(), v2->getLongitude(), v2->getLatitude()) -
                            (Harverstein(v1->getLongitude(), v1->getLatitude(), vset[hubId]->getLongitude(), vset[hubId]->getLatitude()) +
                             Harverstein(v2->getLongitude(), v2->getLatitude(), vset[hubId]->getLongitude(), vset[hubId]->getLatitude()));

            savings.push_back({saving, {v1, v2}});
        }
    }

    // Step 3: Order savings in non-increasing order
    std::sort(savings.rbegin(), savings.rend());

    // Step 4: Apply savings in non-increasing order
    std::unordered_set<Vertex*> visited;
    vset[hubId]->setVisited(true);

    for (const auto& s : savings) {
        Vertex* v1 = s.second.first;
        Vertex* v2 = s.second.second;

        if (v1->isVisited() && v2->isVisited()) // Both already visited, skip
            continue;

        if (!v1->isVisited()) { // v1 visited, add v2 to tour
            tourEdges.push_back(v1->getAdj()[v2->getId()]);
            v2->setVisited(true);
            visited.insert(v2);
        } else if (! v2->isVisited()) {
            tourEdges.push_back(v2->getAdj()[v1->getId()]);
            v1->setVisited(true);
            visited.insert(v1);
        } else {
            Edge* edge1 = v1->getAdj()[hubId];
            Edge* edge2 = v2->getAdj()[hubId];
            if (!edge1 || !edge2)
                continue;

            tourEdges.push_back(edge1);
            tourEdges.push_back(edge2);
            v1->setVisited(true);
            v2->setVisited(true);
            visited.insert(v1);
            visited.insert(v2);
        }
    }

    while (visited.size() > 2) {
        for (auto it = visited.begin(); it != visited.end();) {
            if ((*it)->getId() != hubId) {
                it = visited.erase(it);
            } else {
                ++it;
            }
        }
    }

    return tourEdges;
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

    double cost = INT32_MAX;

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
        for (std::pair<int, Edge*> pair : curr->getAdj()) {
            Edge* e = pair.second;
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

    for(std::pair<int, Edge *> pair: curr->getAdj()){
        Edge* e = pair.second;
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

std::vector<Vertex*> generateRandomTour(Graph& graph) {
    std::vector<Vertex*> tour;

    Vertex* startVertex = nullptr;
    for (const auto& pair : graph.getVertexSet()) {
        if (pair.second->getLatitude() == 0.0 && pair.second->getLongitude() == 0.0) {
            startVertex = pair.second;
            break;
        }
    }

    if (startVertex == nullptr) {
        return tour;
    }

    tour.push_back(startVertex);

    std::unordered_set<Vertex*> visited;
    visited.insert(startVertex);

    while (tour.size() < graph.getVertexSet().size()) {
        bool foundNext = false;
        std::vector<Vertex*> adjacentVertices;

        for(auto pair: startVertex->getAdj()){
            Edge* e = pair.second;
            adjacentVertices.push_back(e->getDestination());
        }

        for (Vertex* vertex : adjacentVertices) {
            if (visited.find(vertex) == visited.end()) {
                tour.push_back(vertex);
                visited.insert(vertex);
                startVertex = vertex;
                foundNext = true;
                break;
            }
        }

        if (!foundNext) {
            break;
        }
    }

    return tour;
}

double evaluateTour(Graph& graph, std::vector<Vertex*>& tour) {
    double cost = 0.0;
    for (size_t i = 0; i < tour.size() - 1; ++i) {
        Vertex* currentVertex = tour[i];
        Vertex* nextVertex = tour[i + 1];
        Edge* edge = currentVertex->getPath();
        if (edge != nullptr) {
            cost += edge->getWeight();
        } else {

            if (currentVertex->getLongitude() == 0 || currentVertex->getLatitude() == 0 ||
                nextVertex->getLongitude() == 0 || nextVertex->getLatitude() == 0) {
                return cost;
            }

            cost += graph.Harverstein(currentVertex->getLongitude(), currentVertex->getLatitude(),nextVertex->getLongitude(), nextVertex->getLatitude());
        }
    }

    Vertex* lastVertex = tour.back();
    Vertex* firstVertex = tour.front();
    Edge* edge = lastVertex->getPath();
    if (edge != nullptr) {
        cost += edge->getWeight();
    } else {
        if (lastVertex->getLongitude() == 0 || lastVertex->getLatitude() == 0 ||
            firstVertex->getLongitude() == 0 || firstVertex->getLatitude() == 0) {
            tour.clear();
            return cost;
        }
        cost += graph.Harverstein(lastVertex->getLongitude(), lastVertex->getLatitude(),firstVertex->getLongitude(), firstVertex->getLatitude());
    }
    return cost;
}




bool kOptMove(Graph& graph, std::vector<Vertex*>& tour, int k,double &curr_cost) {
    bool improvement = false;

    for (size_t i = 0; i < tour.size() - k + 1; ++i) {
        std::vector<Vertex*> newTour = tour;

        newTour.erase(newTour.begin() + i, newTour.begin() + i + k);

        for (size_t j = 0; j < k - 1; ++j) {
            Edge* edge = graph.findEdge(newTour[i + j]->getId(), newTour[i + j + 1]->getId());
            newTour[i + j]->setPath(edge);
        }

        Edge* edge = graph.findEdge(newTour[i + k - 1]->getId(), newTour[(i + k) % tour.size()]->getId());
        newTour[i + k - 1]->setPath(edge);

        double newCost = evaluateTour(graph, newTour);


        if(newTour.empty()){
            improvement = true;
            break;
        }

        if (newCost < evaluateTour(graph, tour)) {
            tour = newTour;
            improvement = true;
            curr_cost = newCost;
        }
    }
    return improvement;
}

std::vector<Vertex*> linKernighan(Graph& graph) {
    std::vector<Vertex*> tour = generateRandomTour(graph);
    double curr_cost = DBL_MAX;

    int max_iterations = 1000;
    int iteration_count = 0;

    while (iteration_count < max_iterations) {
        bool improvement = false;

        for (int k = 2; k <= 3; ++k) {
            improvement |= kOptMove(graph, tour, k, curr_cost);
        }

        if (!improvement) {
            break;
        }

        iteration_count++;
    }

    std::cout << "The cost of the best tour is: " << curr_cost << std::endl;

    return tour;
}
