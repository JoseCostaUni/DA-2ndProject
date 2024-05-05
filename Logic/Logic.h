//Calling different algorithms from here

/**
 * @file Logic.h
 * @brief Declaration of all Logic Functions used during the Project.
 */

#include "../stdafx.h"
#include "Algorithms.h"

/**
 * @brief Retrieves all the source vertices in the given graph.
 * @param g Pointer to the graph.
 */
void getSources(Graph<DeliverySite>* g);

/**
 * @brief Creates a super source and super sink in the graph with connections to all original source and sink nodes.
 * @param g Pointer to the graph.
 * @param SuperSource The super source vertex to be added.
 * @param SuperSink The super sink vertex to be added.
 */
void createSuperSourceSink(Graph<DeliverySite>* g,DeliverySite SuperSource,DeliverySite SuperSink);
/**
 * @brief Removes the super source and super sink from the graph.
 * @param g Pointer to the graph.
 * @param SuperSource The super source vertex to be removed.
 * @param SuperSink The super sink vertex to be removed.
 */
void removeSuperSourceSink(Graph<DeliverySite>* g,DeliverySite SuperSource,DeliverySite SuperSink);
/**
 * @brief Creates a super source in the graph with connections to all original source nodes.
 * @param g Pointer to the graph.
 * @param SuperSource The super source vertex to be added.
 */
void createSuperSource(Graph<DeliverySite>* g,DeliverySite SuperSource);
/**
 * @brief Removes the super source from the graph.
 * @param g Pointer to the graph.
 * @param SuperSource The super source vertex to be removed.
 */
void removeSuperSource(Graph<DeliverySite>* g,DeliverySite SuperSource);

/**
 * @brief Calculates the maximum flow within a single city.
 * @param g Pointer to the graph.
 * @param target The target vertex representing the city.
 */
void calculateMaxFlowInACity(Graph<DeliverySite>* g , DeliverySite& target );

//Find the max flow in the whole graph we use super-source and super-sink
//The super-source node is connected to all the original source nodes by edges with infinite capacity,
// and the super-sink node is connected to all the original sink nodes by edges with capacity equal to their demand

/**
 * @brief Calculates the maximum flow in the entire network using the Ford-Fulkerson algorithm.
 * @param g Pointer to the graph.
 */
void calculateMaxFlowInEntireNetwork(Graph<DeliverySite>* g);

/**
 * @brief Calculates the maximum flow in the entire network using the super source approach.
 * @param g Pointer to the graph.
 * @param target The target vertex representing the city.
 */
void maxFlowWithSuperSource(Graph<DeliverySite>* g , DeliverySite& target);

/**
 * @brief Retrieves all the pipes (edges) in the given graph.
 * @param g Pointer to the graph.
 * @return A vector containing all the edges in the graph.
 */
std::vector<Edge<DeliverySite>*> getPipes(Graph<DeliverySite>* g);