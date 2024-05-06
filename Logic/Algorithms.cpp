#include "Algorithms.h"


std::vector<Vertex*> findTour(Graph * graph , Vertex* startVertex) {
    int numVertices = graph->getNumVertex();

    // Start with any vertex as the initial vertex
    Vertex* currentVertex = startVertex;

    // Initialize a set to keep track of visited vertices
    std::unordered_set<Vertex*, VertexHash, VertexEqual> visited;
    visited.insert(startVertex);

    std::vector<Vertex*> path;
    path.push_back(startVertex);

    // Repeat until all vertices are visited
    for (int i = 0; i < numVertices - 1; ++i) {
        double minDistance = std::numeric_limits<double>::max();
        Vertex* nextVertex = nullptr;

        // Find the nearest unvisited neighbor
        for (Edge* edge : currentVertex->getAdj()) {
            Vertex* neighbor = edge->getDestination();
            if (!visited.count(neighbor) && edge->getWeight() < minDistance) {
                minDistance = edge->getWeight();
                nextVertex = neighbor;
            }
        }

        // Move to the nearest unvisited neighbor
        if (nextVertex != nullptr) {
            visited.insert(nextVertex);
            currentVertex = nextVertex;
            path.push_back(currentVertex);
        }
    }

    // Return to the starting vertex to complete the tour
    std::unordered_set<Edge *> ajdEdjes = currentVertex->getAdj();

    path.push_back(startVertex);

    return path;
}


Vertex* minKey(const Graph* graph ,  const std::vector<double>& key, const std::unordered_set<Vertex*, VertexHash, VertexEqual>& mstSet) {
    double min = std::numeric_limits<double>::max();
    Vertex* minVertex = nullptr;

    for (Vertex* v : graph->getVertexSet()) {
        if (!mstSet.count(v) && key[v->getId()] < min) {
            min = key[v->getId()];
            minVertex = v;
        }
    }

    return minVertex;
}

std::vector<Vertex*> PrimMst(const Graph* graph){
    std::unordered_set<Vertex*, VertexHash, VertexEqual> mstSet;
    std::vector<double> key(graph->getVertexSet().size(), std::numeric_limits<double>::max());
    std::vector<Vertex*> parent(graph->getVertexSet().size(), nullptr);

    // Choose first vertex as the starting point
    key[0] = 0;

    // Loop through all vertices
    for (int count = 0; count < graph->getNumVertex() - 1; ++count) {
        // Get the vertex with the minimum key value from the set of vertices
        Vertex * u = minKey(graph , key, mstSet);

        // Add the selected vertex to the MST set
        mstSet.insert(u);

        // Update key value and parent index of the adjacent vertices of the picked vertex
        for (Edge* edge : u->getAdj()) {
            Vertex* v = edge->getDestination();
            if (!mstSet.count(v) && edge->getWeight() < key[v->getId()]) {
                parent[v->getId()] = u;
                key[v->getId()] = edge->getWeight();
            }
        }
    }
    return parent;
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
   3. Do preorder wack of T. and return Hamilton Cycle
 * */
std::vector<Edge *> TriangularApproximationHeuristic(Graph * graph , Vertex * source , Vertex * dest){

    std::vector<Edge *> optimalRoute;


    return  optimalRoute;
}


