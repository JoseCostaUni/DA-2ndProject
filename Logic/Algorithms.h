#ifndef ALGORITHMS_H
#define ALGORITHMS_H

#include "../Graph.h"
#include <vector>
#include <iostream>
#include <algorithm>
#include <limits>

void tsp(int start, int currPos, int n, double cost, double& minCost, Graph& graph);
double tsp(Graph& graph);

#endif /* ALGORITHMS_H */
