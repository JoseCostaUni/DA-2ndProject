#ifndef ALGORITHMS_H
#define ALGORITHMS_H

/**
 * @file Algorithms.h
 * @brief Declaration of all algorithms used during the Project.
 */

#include "../stdafx.h"


/**
 * @brief Prints a message to the console.
 * @tparam T The type of the message to print.
 * @param _msg The message to print.
 * @param _newline Flag indicating whether to append a newline character.
 */
template <class T>
void print(T _msg , bool _newline){
    std::cout << _msg << " ";
    if(_newline)
        std::cout << std::endl;
}

/**
 * @brief Implements the Edmonds-Karp algorithm to find the maximum flow in the graph.
 * @param g Pointer to the graph.
 * @param source The source delivery site.
 * @param target The target delivery site.
 * @param removed The removed delivery site.
 * @return The maximum flow in the graph.
 */
double edmondsKarp(Graph<DeliverySite> *g, const DeliverySite& source, const DeliverySite& target, const DeliverySite& removed);



double edmondsKarpPipe(Graph<DeliverySite> *g, const DeliverySite& source, const DeliverySite& target,const std::vector<Edge<DeliverySite>*> pump);

/**
 * @brief Calculates the maximum flow through the graph.
 * @param vertexSet The set of vertices in the graph.
 * @return The maximum flow in the graph.
 */
double calculateMaxFlow(std::vector<Vertex<DeliverySite>*>& vertexSet);

/**
 * @brief Calculates the average capacity of pipes in a given set of edges.
 * @param pipes The set of edges representing pipes.
 * @return The average capacity of pipes.
 */
double averagePipeCapacity(const std::vector<Edge<DeliverySite>*>& pipes);

/**
 * @brief Calculates the variance of pipe capacity and flow in a given set of edges.
 * @param pipes The set of edges representing pipes.
 * @param varianceInEachPoint A pointer to a vector to store variance information for each edge.
 * @return The variance of pipe capacity and flow.
 */
double variancePipeCapacityFlow(const std::vector<Edge<DeliverySite>*>& pipes , std::vector<std::pair<double , Edge<DeliverySite>*>>* varianceInEachPoint);

/**
 * @brief Finds the edge with the maximum difference between capacity and flow in a given set of edges.
 * @param pipes The set of edges representing pipes.
 * @return A pair containing the maximum difference and the edge.
 */
std::pair<double , Edge<DeliverySite>*>  maximumDIfferenceCapacityFlow(const std::vector<Edge<DeliverySite>*>& pipes);

/**
 * @brief Calculates the minimum left-over capacity along a given path.
 * @param path The path for which to calculate the minimum left-over capacity.
 * @return The minimum left-over capacity.
 */
double minLeftOverCap(std::vector<Edge<DeliverySite>*>& path);


/**
 * @brief Prints the distance between vertices in the graph.
 * @param g Pointer to the graph.
 */
void printDistance(Graph<DeliverySite>* g);

/**
 * @brief Performs Dijkstra's algorithm to find the shortest paths in the graph.
 * @param g Pointer to the graph.
 * @param root Pointer to the root vertex.
 * @param target Pointer to the target vertex.
 */
void Dijkstra(Graph<DeliverySite>*g , Vertex<DeliverySite>* root , Vertex<DeliverySite>* target);

/**
 * @brief Calculates a path between two vertices in the graph.
 * @param g Pointer to the graph.
 * @param root Pointer to the root vertex.
 * @param target Pointer to the target vertex.
 * @return A pair containing the calculated path and its length.
 */
std::pair<std::vector<Vertex<DeliverySite>*> , int> calculatePath(Graph<DeliverySite>* g , Vertex<DeliverySite>* root , Vertex<DeliverySite>* target);

/**
 * @brief Pumps water along a given path with a specified flow.
 * @param path The path along which to pump water.
 * @param flowToPump The flow rate of water to pump.
 */
void pumpWater(std::vector<Edge<DeliverySite>*>& path , double flowToPump);
//std::vector<Edge<DeliverySite>*> calculatePath(Graph<DeliverySite>* g, Vertex<DeliverySite>* source , Vertex<DeliverySite>* target , Edge<DeliverySite>* edgeToAvoid);

//1º - sort edges por mais flow. se flow igual sort pelas que têm menos capacidade
//2º - encontrar todos os caminhos da origem dessa edge ate ao destino
//3º  - encontrar aquele com maximo leftOver capacity
//4º - dar pump a essa agua
//5º - repeat
/**
 * @brief Applies a heuristic algorithm to optimize pumping stations.
 * @param g Pointer to the graph.
 * @return The metrics after applying the heuristic.
 */
Metrics heuristic(Graph<DeliverySite>*g);

/**
 * @brief Redistributes water without using the Max Flow algorithm.
 *
 * This function redistributes water in the network graph without utilizing the Max Flow algorithm.
 * It iteratively redistributes water from one source to other vertices based on available paths.
 *
 * @param g Pointer to the graph representing the water distribution network.
 * @param removed Pointer to the vertex representing the removed site.
 * @return True if redistribution fails, False otherwise.
 */
bool redistributeWithoutMaxFlowAlgorithm(Graph<DeliverySite>*g, Vertex<DeliverySite>* removed);

/**
 * @brief Finds all possible paths in the graph starting from a given source vertex.
 *
 * This function recursively finds all possible paths starting from the given source vertex.
 *
 * @param g Pointer to the graph representing the water distribution network.
 * @param source Pointer to the source vertex of the paths.
 * @param path Vector storing the current path being explored.
 * @param paths Vector storing all found paths.
 */
void findAllPathsRedistribute(Graph<DeliverySite>*g,Vertex<DeliverySite>* source, std::vector<Edge<DeliverySite>*>& path, std::vector<std::vector<Edge<DeliverySite>*>>& paths);

/**
 * @brief Redistributes water without using the Max Flow algorithm.
 *
 * This function redistributes water in the network graph without utilizing the Max Flow algorithm.
 * It redistributes water along all given paths.
 *
 * @param g Pointer to the graph representing the water distribution network.
 * @param paths Vector of vectors storing paths for water redistribution.
 * @return The total redistributed flow.
 */
double redistributeWaterWithoutMaxFlow2(Graph<DeliverySite>*g, std::vector<std::vector<Edge<DeliverySite>*>>& paths);

/**
 * @brief Creates super source and sink vertices in the graph.
 *
 * This function creates super source and sink vertices and connects them to the appropriate vertices in the graph.
 *
 * @param g Pointer to the graph representing the water distribution network.
 * @param SuperSource The DeliverySite object representing the super source.
 * @param SuperSink The DeliverySite object representing the super sink.
 */
void createSuperSourceSink_(Graph<DeliverySite>* g,DeliverySite SuperSource,DeliverySite SuperSink);
#endif