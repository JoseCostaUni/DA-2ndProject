#include "../Graph.h"

#include <vector>
#include <iostream>
#include <algorithm>
#include <limits>

#define INF std::numeric_limits<double>::max()

void tsp(int start, int currPos, int n, double cost, double& minCost, Graph& graph) {
    Vertex* currVertex = graph.findVertex(currPos);
    Vertex* startVertex = graph.findVertex(start);

    //If all edges have already been visited and there is an edge that leads back to the start vertex
    if (n == 0 && currVertex->getAdj().count(*startVertex->getAdj().begin()) > 0) {
        double cycleCost = cost + (*(startVertex->getAdj().begin()))->getWeight();
        minCost = std::min(minCost, cycleCost);
        return;
    }

    for (auto edge : currVertex->getAdj()) {
        auto nextVertex = edge->getDestination();
        if (!nextVertex->isVisited()) {
            nextVertex->setVisited(true);
            tsp(start, nextVertex->getId(), n - 1, cost + edge->getWeight(), minCost, graph);
            nextVertex->setVisited(false);
        }
    }

}

double tsp(Graph& graph) {
    int n = graph.getNumVertex();
    double minCost = INF;

    for(auto vertex : graph.getVertexSet()){
        vertex->setVisited(false);
    }

    for (auto vertex : graph.getVertexSet()) {
        tsp(vertex->getId(), vertex->getId(), n - 1, 0, minCost, graph);
    }

    return minCost;
}
