#ifndef LOADING_FUNCTIONS_H
#define LOADING_FUNCTIONS_H

/**
 * @file LoadingFunctions.h
 * @brief Declaration of all loading Functions used during the Project.
 */

#include "../stdafx.h"
#include "../Graph.h"

/**
 * @brief Removes termination characters (e.g., '\r', '\n', '\0') from the given string.
 * @param str The string from which to remove termination characters.
 * @complexity O(n), where n is the length of the string.
 */
void Remove_terminations(std::string& str);


double HarvesteinFormula(double latitude1 , double  longitude1 , double latitude2 , double longitude2);

/**
 * @brief Loads cities data from a CSV file.
 * @complexity O(n), where n is the number of lines in the CSV file.
 */
void LoadToyGraphs(Graph * g , const std::string& path , const int& graph);

/**
 * @brief Loads pipes data from a CSV file.
 * @complexity O(n), where n is the number of lines in the CSV file.
 */
void LoadRealWorldGraphs(Graph * g , const std::string& path , const int& graph);

/**
 * @brief Loads water reservoirs data from a CSV file.
 * @complexity O(n), where n is the number of lines in the CSV file.
 */
void LoadMediumGraphs(Graph * g , const std::string& path , const int& graph);

/**
 * @brief Creates a graph based on the loaded data.
 * @param g Pointer to the graph to be created.
 * @return True if the graph creation is successful, false otherwise.
 * @complexity O(m + n), where m is the number of edges (pipes) and n is the number of vertices (delivery sites).
 */
bool createGraph(Graph* g);

#endif // LOADING_FUNCTIONS_H
