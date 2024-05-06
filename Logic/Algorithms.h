
#include "../stdafx.h"

std::vector<Vertex * > NearestNeighbour(Graph * graph , Vertex* startVertex);

Vertex * minKey(const Graph * graph ,  const std::vector<double>& key, const std::unordered_set<Vertex*, VertexHash, VertexEqual>& mstSet);
std::vector<Vertex*> PrimMst(const Graph* graph);

std::vector<Edge *> TriangularApproximationHeuristic(Graph * graph , Vertex * source , Vertex * dest);