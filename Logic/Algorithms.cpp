
/**
 * @file Algorithms.cpp
 * @brief Implementation of all algorithms used during the Project.
 */

#include <climits>
#include "Algorithms.h"
#include "../stdafx.h"
#include "Logic.h"



/**
 * @brief Finds the minimum residual capacity along an augmenting path from a source vertex to a sink vertex.
 * @param source The source vertex.
 * @param sink The sink vertex.
 * @return The minimum residual capacity along the augmenting path.
 * @details Time Complexity: O(N), where N is the number of vertices along the path.
 */
double findMinResidualAlongPath(Vertex<DeliverySite> *source, Vertex<DeliverySite> *sink) {
    double f = DBL_MAX; // Traverse the augmenting path to find the minimum residual capacity

    for (Vertex<DeliverySite>* v = sink; v != source;) {
        Edge<DeliverySite>* e = v->getPath();
        if (e->getDest() == v){
            f = std::min(f, e->getWeight() - e->getFlow());
            v = e->getOrig();
        } else {
            f = std::min(f, e->getFlow());
            v = e->getDest();
        }
    }
    return f; // Return the minimum residual capacity
}

/**
 * @brief Finds an augmenting path in the graph using BFS.
 * @param g Pointer to the graph.
 * @param source Pointer to the source vertex.
 * @param sink Pointer to the sink vertex.
 * @param removed Pointer to the removed vertex (optional).
 * @return True if an augmenting path is found, otherwise false.
 * @details Time Complexity: O(V + E), where:
 * - V is the number of vertices in the graph.
 * - E is the number of edges in the graph.
 */
bool findAugmentingPath(Graph<DeliverySite> *g, Vertex<DeliverySite> *source, Vertex<DeliverySite> *sink, Vertex<DeliverySite> *removed) {
// Mark all vertices as not visited
    for(Vertex<DeliverySite>* v : g->getVertexSet()) {
        v->setVisited(false);
    }
// Mark the source vertex as visited and enqueue it

    if (removed != nullptr) removed->setVisited(true);

    std::queue<Vertex<DeliverySite> *> q;
    q.push(source);
    source->setVisited(true);

// BFS to find an augmenting path
    while(!q.empty() && !sink->isVisited()) {
        Vertex<DeliverySite>* v = q.front();
        q.pop();
        for (Edge<DeliverySite>* e : v->getAdj()) {
            if(e->isSelected()){
                continue;
            }
            Vertex<DeliverySite>* dest = e->getDest();
            if (!dest->isVisited() && (e->getWeight() - e->getFlow() > 0)) {
                dest->setVisited(true);
                dest->setPath(e);
                q.push(dest);
            }
        }
        for (Edge<DeliverySite>* e: v->getIncoming()) {
            Vertex<DeliverySite>* origin = e->getOrig();
            if (!origin->isVisited() && (e->getFlow() > 0)) {
                origin->setVisited(true);
                origin->setPath(e);
                q.push(origin);
            }
        }
    }
    return sink->isVisited();
}

/**
 * @brief Augments the flow along the augmenting path found using BFS.
 * @param source Pointer to the source vertex.
 * @param sink Pointer to the sink vertex.
 * @param f The flow to augment along the path.
 * @details Time Complexity: O(N), where N is the number of vertices along the path.
 */
void augmentFlowAlongPath(Vertex<DeliverySite> *source, Vertex<DeliverySite> *sink, double f) {
// Traverse the augmenting path and update the flow values accordingly

    for (Vertex<DeliverySite>* v = sink; v != source;) {
        Edge<DeliverySite>* e = v->getPath();

        double flow = e->getFlow();
        if (e->getDest() == v) {
            e->setFlow(flow + f);
            v = e->getOrig();
        }
        else {
            e->setFlow(flow - f);
            v = e->getDest();
        }
    }


}

/**
 * @brief Implements the Edmonds-Karp algorithm to find the maximum flow in the graph.
 * @param g Pointer to the graph.
 * @param source The source delivery site.
 * @param target The target delivery site.
 * @param removed The removed delivery site.
 * @return The maximum flow in the graph.
 * @details Time Complexity: O((V^2) * E), where:
 * - V is the number of vertices in the graph.
 * - E is the number of edges in the graph.
 */
double edmondsKarp(Graph<DeliverySite> *g, const DeliverySite& source, const DeliverySite& target,const DeliverySite& removed) {
    double maxFlow = 0;
// Find source and target vertices in the graph
    Vertex<DeliverySite>* s = g->findVertex(source);
    Vertex<DeliverySite>* t = g->findVertex(target);
    Vertex<DeliverySite>* remove = g->findVertex(removed);
// Validate source and target vertices
    if (s == nullptr || t == nullptr || s == t)
        throw std::logic_error("Invalid source and/or target vertex");
// Initialize flow on all edges to 0
    for (auto v : g->getVertexSet()) {
        for (auto e: v->getAdj()) {
            e->setFlow(0);
            e->setSelected(false);
        }
    }

// While there is an augmenting path, augment the flow along the path
    while(findAugmentingPath(g, s, t, remove) ) {
        double f = findMinResidualAlongPath(s, t);
        maxFlow += f;
        augmentFlowAlongPath(s, t, f);
    }
    return maxFlow;
}

/**
 * @brief Implements the Edmonds-Karp algorithm to find the maximum flow in the graph without a certain edge.
 * @param g Pointer to the graph.
 * @param source The source delivery site.
 * @param target The target delivery site.
 * @param removed The removed egde.
 * @return The maximum flow in the graph.
 * @details Time Complexity: O((V^2) * E), where:
 * - V is the number of vertices in the graph.
 * - E is the number of edges in the graph.
 */
double edmondsKarpPipe(Graph<DeliverySite> *g, const DeliverySite& source, const DeliverySite& target,const std::vector<Edge<DeliverySite>*> pump) {
    double maxFlow = 0;
// Find source and target vertices in the graph
    Vertex<DeliverySite>* s = g->findVertex(source);
    Vertex<DeliverySite>* t = g->findVertex(target);
// Validate source and target vertices
    if (s == nullptr || t == nullptr || s == t)
        throw std::logic_error("Invalid source and/or target vertex");
// Initialize flow on all edges to 0
    for (auto v : g->getVertexSet()) {
        for (auto e: v->getAdj()) {
            e->setFlow(0);
            e->setSelected(false);
            for(auto edge : pump){
                if(e->getOrig()->getInfo().getCode() == edge->getOrig()->getInfo().getCode() && e->getDest()->getInfo().getCode() == edge->getDest()->getInfo().getCode()){
                    e->setSelected(true);
                }
            }
        }
    }

// While there is an augmenting path, augment the flow along the path
    while(findAugmentingPath(g, s, t, nullptr) ) {
        double f = findMinResidualAlongPath(s, t);
        maxFlow += f;
        augmentFlowAlongPath(s, t, f);
    }
    return maxFlow;
}

//root will be the origin of edge with the lowest flow
//this edge will be locked to guarantee that it is not picked during this algorithm
//we can try to tell this algo to
//djikstra picks the paths with full edges
/**
 * @brief Performs Dijkstra's algorithm to find the shortest paths in terms of flow from a root vertex to a target vertex.
 * @param g Pointer to the graph.
 * @param root Pointer to the root vertex.
 * @param target Pointer to the target vertex.
 * @details Time Complexity: O((V + E) * log(V)), where:
 * - V is the number of vertices in the graph.
 * - E is the number of edges in the graph.
 */
void Dijkstra(Graph<DeliverySite>*g , Vertex<DeliverySite>* root , Vertex<DeliverySite>* target) {

    MutablePriorityQueue<Vertex<DeliverySite>> vertexQueue;
    for(Vertex<DeliverySite>* v : g->getVertexSet()){
        vertexQueue.insert(v);
    }

    for (Vertex<DeliverySite> *v: g->getVertexSet()) {
        v->setDist(INF);
        v->setPath(nullptr);
    }

    root->setDist(0);

    while (!vertexQueue.empty()){
        Vertex<DeliverySite>* u = vertexQueue.extractMin();

        for(Edge<DeliverySite>* e : u->getAdj()){
            //we want the minimum distance by flow
            Vertex<DeliverySite>* v = e->getDest();
            if(v->getDist() > (u->getDist() + e->getFlow())){
                v->setDist(u->getDist() + e->getFlow());
                v->setPath(e);
            }
        }
    }

    target->setDist(INF);
}

/**
 * @brief Calculates the minimum left-over capacity along a given path.
 * @param path The path for which to calculate the minimum left-over capacity.
 * @return The minimum left-over capacity.
 * @details Time Complexity: O(N), where N is the number of edges in the path.
 */
double minLeftOverCap(std::vector<Edge<DeliverySite>*>& path){
    auto a = 0;
    for(auto e : path){
        if(e->getFlow() != e->getWeight())
            a++;
    }

    auto min = DBL_MAX;
    for(Edge<DeliverySite>* edge : path){
        if(edge->getWeight() == edge->getFlow())
            return 0;
        if(edge->getWeight() - edge->getFlow() < min){
            min = edge->getWeight() - edge->getFlow();
        }
    }

    return min;
}

/**
 * @brief Computes metrics heuristic for optimizing water delivery.
 * @param g Pointer to the graph of delivery sites.
 * @return Metrics after optimization.
 * @details Time Complexity: O(((V+E) * E +E) * P), where:
 * - N is the number of delivery sites in the graph.
 * - E is the number of edges in the graph.
 * - P is the number of times that the while cycle is repeated
 */
Metrics heuristic(Graph<DeliverySite>*g){
    std::vector<Edge<DeliverySite>*> edges;

    edges = g->getEdges(); //O(V+E)

    Metrics finalMetrics = g->calculateMetrics();
    Metrics initialMetrics = finalMetrics;

    g->printMetrics(initialMetrics);
    initialMetrics = {DBL_MAX , DBL_MAX , DBL_MAX , DBL_MAX};
    while(finalMetrics.variance < initialMetrics.variance || finalMetrics.avg < initialMetrics.avg){

        std::sort(edges.begin(), edges.end(), [](Edge<DeliverySite>* a, Edge<DeliverySite>* b) {

            if(a->getWeight() - a->getFlow() == b->getWeight() - b->getFlow()){
                return a->getWeight() > b->getWeight();
            }

            return a->getWeight() - a->getFlow() < b->getWeight() - b->getFlow();
        });//O(E log E)

        //O(E)
        for(Edge<DeliverySite>* e : edges){
            std::vector<Edge<DeliverySite>*> path;
            std::vector<std::vector<Edge<DeliverySite>*>> allPaths;

            //O(V+E)
            allPaths = g->allPaths(e->getOrig()->getInfo() , e->getDest()->getInfo());

            double maxDiff = -1;
            if(allPaths.empty())
                continue;

            for(std::vector<Edge<DeliverySite>*> tempPath : allPaths){
                double minFlow = minLeftOverCap(tempPath);
                if (minFlow > maxDiff) {
                    maxDiff = minFlow;
                    path = tempPath;
                }
            }

            double waterToPump = maxDiff;
            if(e->getFlow() - waterToPump < 0)
                waterToPump = e->getFlow();

            e->setFlow(e->getFlow() - waterToPump);

            pumpWater(path , waterToPump);
        }

        initialMetrics = finalMetrics;

        finalMetrics = g->calculateMetrics();

    }

    finalMetrics = g->calculateMetrics();

    g->printMetrics(finalMetrics);

    for(auto e : g->getEdges()){
        if(e->getFlow() > e->getWeight()){
            print("SOBRECARGAAAA" , false);
        }
        if(e->getFlow() < 0)
            print("DESCEU O CANOOOOO" , false);
    }
    return  finalMetrics;
}

/**
 * @brief Pumps water along the given path with the specified flow rate.
 * @param path The path along which water is to be pumped.
 * @param flowToPump The flow rate of water to pump.
 * @details Time Complexity: O(N), where N is the number of edges in the path.
 */
void pumpWater(std::vector<Edge<DeliverySite>*>& path , double flowToPump){
    if(flowToPump != 0)
        auto a = 0;
    for(Edge<DeliverySite>* e : path){
            e->setFlow(e->getFlow() + flowToPump);
    }
}

/**
 * @brief Redistributes water without using the Max Flow algorithm.
 *
 * This function redistributes water in the network graph without utilizing the Max Flow algorithm.
 * It iteratively redistributes water from one source to other vertices based on available paths.
 *
 * @param g Pointer to the graph representing the water distribution network.
 * @param removed Pointer to the vertex representing the removed site.
 * @return True if redistribution fails, False otherwise.
 *
 * @note This function modifies the flow in the edges of the graph.
 *
 * @see findAugPath, minResidualAugPath, augmentFlowPath
 *
 * @timecomplexity O(E * V^2), where E is the number of edges and V is the number of vertices in the graph.
 */
bool redistributeWithoutMaxFlowAlgorithm(Graph<DeliverySite>*g, Vertex<DeliverySite>* removed){
    bool flag_edmonds_Karp = false;
    std::vector<Vertex<DeliverySite>*> cities;

    std::stack<Vertex<DeliverySite>*> pumpUsed;

    std::stack<std::vector<std::vector<Edge<DeliverySite>*>>> allPaths;

    for(Vertex<DeliverySite>* ver: g->getVertexSet()){
        if(ver->getInfo().getNodeType() == CITY){
            cities.push_back(ver);
        }
    }

    for(Vertex<DeliverySite>* city : cities){
        allPaths.push(g->allPaths(removed->getInfo(),city->getInfo()));
    }

    while(!allPaths.empty()){
        std::vector<std::vector<Edge<DeliverySite>*>> stackTop = allPaths.top();
        allPaths.pop();

        for(std::vector<Edge<DeliverySite>*> path : stackTop){
            double minFlow = DBL_MAX;
            for(auto & i : path) minFlow = std::min(minFlow,i->getFlow());

            for(auto & i : path) i->setFlow(i->getFlow() - minFlow);
        }

    }

    for(Edge<DeliverySite>* edge : g->getEdges()){
        edge->setSelected(false);
    }

    for(Edge<DeliverySite>* edge : removed->getAdj()){
        edge->getDest()->getInfo().setDemand(edge->getFlow());
        pumpUsed.push(edge->getDest());
        edge->setSelected(true);
    }

    while(!pumpUsed.empty()){
        Vertex<DeliverySite>* pump = pumpUsed.top();
        pumpUsed.pop();

        for(Vertex<DeliverySite>* city : cities){
            while(findAugmentingPath(g, pump, city, removed) ) {
                double f = findMinResidualAlongPath(pump, city);
                if(pump->getInfo().getDemand() - f >= 0){
                    augmentFlowAlongPath(pump, city, f);
                    pump->getInfo().setDemand(pump->getInfo().getDemand() - f);
                }else{
                    augmentFlowAlongPath(pump,city, f - abs(pump->getInfo().getDemand() - f));
                    pump->getInfo().setDemand(0);
                    break;
                }
            }

        }

        if(pump->getInfo().getDemand() > 0){
            flag_edmonds_Karp = true;
            break;
        }

    }

    if(flag_edmonds_Karp){
        for(Vertex<DeliverySite>* ver: g->getVertexSet()){
            if(ver->getInfo().getNodeType() == FIRE_STATION){
                ver->getInfo().setDemand(0);
            }
        }
        std::cout << "Failed to find the whole impact" << std::endl;
    }else{
        std::cout << "Impact well found" << std::endl;
    }


    return flag_edmonds_Karp;
}

/**
 * @brief Finds an augmenting path in the graph using BFS.
 *
 * This function finds an augmenting path from the given source vertex to the removed site.
 * It utilizes BFS (Breadth-First Search) algorithm to traverse the graph.
 *
 * @param g Pointer to the graph representing the water distribution network.
 * @param source Pointer to the source vertex.
 * @param removed Pointer to the vertex representing the removed site.
 * @return Pointer to the vertex representing the augmenting path, or nullptr if not found.
 *
 * @timecomplexity O(V + E), where V is the number of vertices and E is the number of edges in the graph.
 */
Vertex<DeliverySite>* findAugPath(Graph<DeliverySite>*g, Vertex<DeliverySite> *source, Vertex<DeliverySite>* removed){
    if(source->getInfo().getNodeType() == CITY && source->getInNeed() == 0) return nullptr;

    for(Vertex<DeliverySite>* v : g->getVertexSet()){
        v->setVisited(false);
    }

    removed->setVisited(true);

    std::queue<Vertex<DeliverySite>*> q;
    q.push(source);
    source->setVisited(true);

    while(!q.empty()) {
        Vertex<DeliverySite>* v = q.front();
        q.pop();
        if(v->getInfo().getNodeType() == WATER_RESERVOIR) {
            return v;
        }

        for(Edge<DeliverySite>* edge : v->getIncoming()){
            Vertex<DeliverySite>* origin = edge->getOrig();
            if(!origin->isVisited() && (edge->getWeight() - edge->getFlow() > 0)){
                if(origin->getInfo().getNodeType() == WATER_RESERVOIR){
                    double remain = origin->getInfo().getMaxDelivery() - origin->calculateOutgoingFlow();
                    if(remain <= 0){
                        origin->setVisited(true);
                        continue;
                    }
                }
                origin->setVisited(true);
                origin->setPath(edge);
                q.push(origin);
            }
        }
    }
    return nullptr;
}

/**
 * @brief Finds the minimum residual along an augmenting path.
 *
 * This function calculates the minimum residual capacity along an augmenting path.
 *
 * @param g Pointer to the graph representing the water distribution network.
 * @param source Pointer to the source vertex of the path.
 * @param sink Pointer to the sink vertex of the path.
 * @return The minimum residual capacity along the augmenting path.
 *
 * @timecomplexity O(V), where V is the number of vertices in the augmenting path.
 */
double minResidualAugPath(Graph<DeliverySite>*g,Vertex<DeliverySite>* source, Vertex<DeliverySite>* sink){
    double maxFlow = DBL_MAX;
    for(Vertex<DeliverySite>* v = source; v != sink;){
        Edge<DeliverySite>* e = v->getPath();

        maxFlow = std::min(maxFlow,e->getWeight() - e->getFlow());

        if(v->getInfo().getNodeType() == WATER_RESERVOIR){
            double remain = v->getInfo().getMaxDelivery() - v->calculateOutgoingFlow();
            maxFlow = std::min(maxFlow,remain);
        }
        v = e->getDest();
    }
    return maxFlow;
}

/**
 * @brief Augments flow along an augmenting path.
 *
 * This function augments the flow along the given augmenting path in the graph.
 *
 * @param source Pointer to the source vertex of the path.
 * @param sink Pointer to the sink vertex of the path.
 * @param flow The amount of flow to be augmented along the path.
 *
 * @timecomplexity O(V), where V is the number of vertices in the augmenting path.
 */
void augmentFlowPath(Vertex<DeliverySite>* source,Vertex<DeliverySite>* sink, double flow){
    for(Vertex<DeliverySite>* v = source; v != sink;){
        Edge<DeliverySite>* edge = v->getPath();

        if(edge->getNeeds() != DBL_MAX) edge->setNeeds(edge->getNeeds() - flow);
        edge->setFlow(edge->getFlow() + flow);
        v = edge->getDest();
    }
}


/**
 * @brief Finds all possible paths in the graph starting from a given source vertex.
 *
 * This function recursively finds all possible paths starting from the given source vertex.
 *
 * @param g Pointer to the graph representing the water distribution network.
 * @param source Pointer to the source vertex of the paths.
 * @param path Vector storing the current path being explored.
 * @param paths Vector storing all found paths.
 *
 * @timecomplexity Exponential in the size of the graph.
 */
void findAllPathsRedistribute(Graph<DeliverySite>*g,Vertex<DeliverySite>* source, std::vector<Edge<DeliverySite>*>& path, std::vector<std::vector<Edge<DeliverySite>*>>& paths){
    source->setVisited(true);

    if(source->getInfo().getNodeType() == CITY){
        paths.push_back(path);
    }else{
        for(Edge<DeliverySite>* edge : source->getAdj()){
            if(!edge->getDest()->isVisited()){
                path.push_back(edge);
                findAllPathsRedistribute(g,edge->getDest(),path,paths);
            }
        }
    }

    if(!path.empty()) path.pop_back();
    source->setVisited(false);
}


/**
 * @brief Redistributes water without using the Max Flow algorithm.
 *
 * This function redistributes water in the network graph without utilizing the Max Flow algorithm.
 * It redistributes water along all given paths.
 *
 * @param g Pointer to the graph representing the water distribution network.
 * @param paths Vector of vectors storing paths for water redistribution.
 * @return The total redistributed flow.
 *
 * @timecomplexity O(E * V^2), where E is the number of edges and V is the number of vertices in the graph.
 */
double redistributeWaterWithoutMaxFlow2(Graph<DeliverySite>*g, std::vector<std::vector<Edge<DeliverySite>*>>& paths){
    if(paths.empty()) return 0;

    for(auto v : g->getVertexSet()){
        v->setInNeed(0);
        v->setAlreadyHas(0);
        v->setVisited(false);
        for(Edge<DeliverySite>* edge: v->getAdj()){
            edge->setSelected(false);
            edge->setNeeds(0);
        }
    }

    std::vector<Vertex<DeliverySite>*> cities;
    for(auto path : paths){
        if(path.size() < 2) continue;
        Edge<DeliverySite>* edge = path[path.size() - 1];
        if(!edge->isSelected()){
            Vertex<DeliverySite>* cityFound = edge->getDest();
            cityFound->setInNeed(cityFound->getInNeed() + edge->getFlow());
            edge->setNeeds(edge->getFlow());
            cities.push_back(edge->getDest());
        }
    }

    for(std::vector<Edge<DeliverySite>*> path : paths){
        double minFlow = DBL_MAX;
        for(auto & i : path) minFlow = std::min(minFlow,i->getFlow());

        for(auto & i : path) i->setFlow(i->getFlow() - minFlow);
    }

    for(Vertex<DeliverySite>* v: g->getVertexSet()){
        for(Edge<DeliverySite>* e : v->getAdj()){
            e->setSelected(false);
        }
    }

    for(std::vector<Edge<DeliverySite>*> path : paths){
        if(path.size() < 2) continue;
        Edge<DeliverySite>* edge = path[path.size() - 1];
        if(!edge->isSelected()){
            Vertex<DeliverySite>* city = edge->getDest();
            city->setInNeed(city->getInNeed() - edge->getFlow());
            edge->setNeeds(edge->getNeeds() - edge->getFlow());
            edge->setSelected(true);
        }
    }

    Vertex<DeliverySite>* water_reservoir = paths[0][0]->getOrig();

    for(Vertex<DeliverySite>* city : cities){
        while(true){
            Vertex<DeliverySite>* augment = findAugPath(g,city,water_reservoir);
            if(augment == nullptr){
                break;
            }
            double flow = minResidualAugPath(g,augment,city);

            if(flow > city->getInNeed()){
                flow = city->getInNeed();
            }
            augmentFlowPath(augment,city,flow);
            city->setInNeed(city->getInNeed() - flow);
            city->setAlreadyHas(city->getAlreadyHas() + flow);

        }
    }

    double flow = 0;
    for(Vertex<DeliverySite>* v : g->getVertexSet()){
        if(v->getInfo().getNodeType() == CITY ){
            for(auto e : v->getIncoming()){
                flow += e->getFlow();
            }
        }
    }

    return flow;
}


/**
 * @brief Creates super source and sink vertices in the graph.
 *
 * This function creates super source and sink vertices and connects them to the appropriate vertices in the graph.
 *
 * @param g Pointer to the graph representing the water distribution network.
 * @param SuperSource The DeliverySite object representing the super source.
 * @param SuperSink The DeliverySite object representing the super sink.
 *
 * @timecomplexity O(V + E), where V is the number of vertices and E is the number of edges in the graph.
 */
void createSuperSourceSink_(Graph<DeliverySite>* g,DeliverySite SuperSource,DeliverySite SuperSink){
    g->addVertex(SuperSource);
    g->addVertex(SuperSink);

    for (auto v: g->getVertexSet()) {
        nodeTypes code = v->getInfo().getNodeType();
        if (code == WATER_RESERVOIR){
            g->addEdge(SuperSource,v->getInfo(),v->getInfo().getMaxDelivery());
        } else if (code == CITY && v->getInfo().getCode() != "SuperSource" && v->getInfo().getCode() != "SuperSink") {
            g->addEdge(v->getInfo(),SuperSource, v->getInfo().getDemand());
        }
    }

    for (auto v : g->getVertexSet()) {
        if (v->getInfo().getCode() != "SuperSource" && v->getInfo().getCode() != "SuperSink" && v->getInfo().getNodeType() == CITY) {
            g->addEdge(v->getInfo(), SuperSink, v->getInfo().getDemand());
        }
    }

    for(auto v : g->getVertexSet()){
        v->setVisited(false);
        for(Edge<DeliverySite>* edge: v->getAdj()){
            edge->setSelected(false);
        }
    }

}
