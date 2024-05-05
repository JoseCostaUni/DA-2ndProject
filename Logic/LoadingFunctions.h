#ifndef LOADING_FUNCTIONS_H
#define LOADING_FUNCTIONS_H

/**
 * @file LoadingFunctions.h
 * @brief Declaration of all loading Functions used during the Project.
 */

#include <unordered_set>
#include <vector>
#include "DeliverySites.h"
#include "PumpingStations.h"
#include "../data_structures/Graph.h"

/**
 * @brief Normalizes two strings by removing leading and trailing whitespaces and converting them to lowercase.
 * @param str1 The first string to be normalized.
 * @param str2 The second string to be normalized.
 * @complexity O(n), where n is the length of the longer string between str1 and str2.
 */
void NormalizeString(std::string& str1 , std::string& str2);

/**
 * @brief Removes termination characters (e.g., '\r', '\n', '\0') from the given string.
 * @param str The string from which to remove termination characters.
 * @complexity O(n), where n is the length of the string.
 */
void Remove_terminations(std::string& str);

/**
 * @brief Loads cities data from a CSV file.
 * @complexity O(n), where n is the number of lines in the CSV file.
 */
void LoadCities(const std::string& path);

/**
 * @brief Loads pipes data from a CSV file.
 * @complexity O(n), where n is the number of lines in the CSV file.
 */
void LoadPipes(const std::string& path);

/**
 * @brief Loads water reservoirs data from a CSV file.
 * @complexity O(n), where n is the number of lines in the CSV file.
 */
void LoadWaterReservoirs(const std::string& path);

/**
 * @brief Loads fire stations data from a CSV file.
 * @complexity O(n), where n is the number of lines in the CSV file.
 */
void LoadFireStations(const std::string& path);

/**
 * @brief Creates a graph based on the loaded data.
 * @param g Pointer to the graph to be created.
 * @return True if the graph creation is successful, false otherwise.
 * @complexity O(m + n), where m is the number of edges (pipes) and n is the number of vertices (delivery sites).
 */
bool createGraph(Graph<DeliverySite>* g);

#endif // LOADING_FUNCTIONS_H
