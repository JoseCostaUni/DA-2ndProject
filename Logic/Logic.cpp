#include <climits>
#include "Logic.h"

/**
 * @file Logic.cpp
 * @brief Implementation of all Logic Functions used during the Project.
 */

std::vector<Vertex<DeliverySite>*> sources;
std::vector<Vertex<DeliverySite>*> sinks;
std::unordered_map<Vertex<DeliverySite>* , int> flowMap;

/**
 * @brief Retrieves all the source nodes (water reservoirs) from the graph.
 * @param g Pointer to the graph containing the delivery sites.
 * @complexity O(n), where n is the number of vertices in the graph.
 */
void getSources(Graph<DeliverySite>* g){
    for(Vertex<DeliverySite>* deliverySite : g->getVertexSet()){
        if(deliverySite->getInfo().getNodeType() == WATER_RESERVOIR){
            sources.push_back(deliverySite);
        }
    }
}

/**
 * @brief Retrieves all the sink nodes (cities) from the graph.
 * @param g Pointer to the graph containing the delivery sites.
 * @complexity O(n), where n is the number of vertices in the graph.
 */
void getSinks(Graph<DeliverySite>* g){
    for(Vertex<DeliverySite>* deliverySite : g->getVertexSet()){
        if(deliverySite->getInfo().getNodeType() == CITY){
            sinks.push_back(deliverySite);
        }
    }
}

//Find the max flow in the whole graph we use super-source and super-sink
//The super-source node is connected to all the original source nodes by edges with infinite capacity,
// and the super-sink node is connected to all the original sink nodes by edges with capacity equal to their demand
/**
 * @brief Creates a super-source and super-sink for the graph, connecting them to appropriate nodes.
 * The super-source is connected to all the original source nodes with infinite capacity,
 * and the super-sink is connected to all the original sink nodes with capacity equal to their demand.
 * @param g Pointer to the graph containing the delivery sites.
 * @param SuperSource DeliverySite object representing the super-source node.
 * @param SuperSink DeliverySite object representing the super-sink node.
 * @complexity O(n), where n is the number of vertices in the graph.
 */
void createSuperSourceSink(Graph<DeliverySite>* g,DeliverySite SuperSource,DeliverySite SuperSink){
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

}

/**
 * @brief Creates a super-source for the graph, connecting it to appropriate nodes.
 * The super-source is connected to all the original source nodes with infinite capacity.
 * @param g Pointer to the graph containing the delivery sites.
 * @param SuperSource DeliverySite object representing the super-source node.
 * @complexity O(n), where n is the number of vertices in the graph.
 */
void createSuperSource(Graph<DeliverySite>* g,DeliverySite SuperSource){
    g->addVertex(SuperSource);

    for (auto v: g->getVertexSet()) {
        nodeTypes code = v->getInfo().getNodeType();
        if (code == WATER_RESERVOIR){
            g->addEdge(SuperSource,v->getInfo(),v->getInfo().getMaxDelivery());
        } else if (code == CITY && v->getInfo().getCode() != "SuperSource" && v->getInfo().getCode() != "SuperSink") {
            g->addEdge(v->getInfo(),SuperSource, v->getInfo().getDemand());
        }
    }

}

/**
 * @brief Removes the super-source and super-sink nodes from the graph.
 * @param g Pointer to the graph containing the delivery sites.
 * @param SuperSource DeliverySite object representing the super-source node.
 * @param SuperSink DeliverySite object representing the super-sink node.
 * @complexity O(n), where n is the number of vertices in the graph.
 */
void removeSuperSourceSink(Graph<DeliverySite>* g,DeliverySite SuperSource,DeliverySite SuperSink) {
    g->removeVertex(SuperSource);
    g->removeVertex(SuperSink);
}

/**
 * @brief Removes the super-source node from the graph.
 * @param g Pointer to the graph containing the delivery sites.
 * @param SuperSource DeliverySite object representing the super-source node.
 * @complexity O(n), where n is the number of vertices in the graph.
 */
void removeSuperSource(Graph<DeliverySite>* g,DeliverySite SuperSource) {
    g->removeVertex(SuperSource);
}

/**
 * @brief Retrieves all the unique pipes (edges) from the graph.
 * @param g Pointer to the graph containing the delivery sites.
 * @return A vector of pointers to unique edges in the graph.
 * @complexity O(m), where m is the number of edges in the graph.
 */
std::vector<Edge<DeliverySite>*> getPipes(Graph<DeliverySite>* g){
    std::unordered_set<Edge<DeliverySite>*> res;

    for(auto v : g->getVertexSet()){
        for(Edge<DeliverySite>* e : v->getAdj()){
            if(res.find(e) == res.end())
                res.insert(e);
        }
    }

    std::vector<Edge<DeliverySite>*> result(res.begin(), res.end());

    return result;
}