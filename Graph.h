// Original code by Gonçalo Leão
// Updated by DA 2023/2024 Team

#ifndef DA_TP_CLASSES_GRAPH
#define DA_TP_CLASSES_GRAPH

#include <iostream>
#include <vector>
#include <queue>
#include <limits>
#include <algorithm>
#include <cfloat>
#include <cstdint>
#include <unordered_set>


struct Metrics{
    double avg;
    double variance;
    double maxDiff;
};


class Edge;

#define INF std::numeric_limits<double>::max()

/************************* Vertex  **************************/


class Vertex {
public:

    Vertex(long id, double longitude, double latitude);
    Vertex(long id);
    long getId() const;
    double getLongitude() const;
    double getLatitude() const;
    std::unordered_set<Edge*> getAdj() const;
    std::unordered_set<Edge*> getInc() const;
    void addAdjEdge(Edge* edge);
    void addIncEdge(Edge* edge);
    void deleteIncEdge(Edge* edge);
    void deleteAdjEdge(Edge* edge);

    void setId(long newId);
    void setLongitude(double newLongitude);
    void setLatitude(double newLatitude);
protected:
    long id;
    double longitude;
    double latitude;
    std::unordered_set<Edge*> incomingEdges;
    std::unordered_set<Edge*> adjacentEdges;
    bool visited = false;
    bool processing = false;
    unsigned int indegree = 0;
    double dist = 0;

};

inline Vertex::Vertex(long id, double longitude, double latitude)
        : id(id), longitude(longitude), latitude(latitude) {}

inline Vertex::Vertex(long id) : id(id) {}

inline long Vertex::getId() const {
    return id;
}

inline double Vertex::getLongitude() const {
    return longitude;
}

inline double Vertex::getLatitude() const {
    return latitude;
}

inline std::unordered_set<Edge*> Vertex::getAdj() const {
    return adjacentEdges;
}

inline std::unordered_set<Edge*> Vertex::getInc() const {
    return incomingEdges;
}

inline void Vertex::addAdjEdge(Edge* edge) {
    adjacentEdges.insert(edge);
}

inline void Vertex::addIncEdge(Edge* edge) {
    incomingEdges.insert(edge);
}

inline void Vertex::deleteIncEdge(Edge* edge) {
    incomingEdges.erase(edge);
}

inline void Vertex::deleteAdjEdge(Edge* edge) {
    adjacentEdges.erase(edge);
}

inline void Vertex::setId(long newId) {
    id = newId;
}

inline void Vertex::setLongitude(double newLongitude) {
    longitude = newLongitude;
}

inline void Vertex::setLatitude(double newLatitude) {
    latitude = newLatitude;
}


struct VertexHash {
    std::size_t operator()(const Vertex* v) const {
        // Hash the vertex ID
        return std::hash<long>{}(v->getId());
    }
};

// Custom equality operator for Vertex
struct VertexEqual {
    bool operator()(const Vertex* lhs, const Vertex* rhs) const {
        // Compare vertex IDs for equality
        return lhs->getId() == rhs->getId();
    }
};

/********************** Edge  ****************************/


class Edge {
public:
    Edge(Vertex* source, Vertex* destination, double weight);
    Vertex* getSource() const;
    Vertex* getDestination() const;
    double getWeight() const;

    void setSource(Vertex* newSource);
    void setDestination(Vertex* newDestination);
    void setWeight(double newWeight);


private:
    Vertex* source;
    Vertex* destination;
    double weight;
};

/********************** Graph  ****************************/


class Graph {
public:
    ~Graph();
    Vertex* findVertex(long id) const;
    bool addVertex(const Vertex* vertex);
    bool removeVertex(long id);
    bool addEdge(long sourceId, long destId, double weight);
    bool removeEdge(long sourceId, long destId);
    bool addBidirectionalEdge(long sourceId, long destId, double weight);
    int getNumVertex() const;
    std::unordered_set<Vertex* , VertexHash, VertexEqual> getVertexSet() const;

    void clear();
    void printNodesContente() const;
private:
    std::unordered_set<Edge*> edgeSet;
    std::unordered_set<Vertex* , VertexHash, VertexEqual> vertexSet;
};

/********************** Edge  ****************************/


inline Edge::Edge(Vertex* source, Vertex* destination, double weight)
        : source(source), destination(destination), weight(weight) {}

inline Vertex* Edge::getSource() const {
    return source;
}

inline Vertex* Edge::getDestination() const {
    return destination;
}

inline double Edge::getWeight() const {
    return weight;
}

inline void Edge::setSource(Vertex* newSource) {
    source = newSource;
}

inline void Edge::setDestination(Vertex* newDestination) {
    destination = newDestination;
}

inline void Edge::setWeight(double newWeight) {
    weight = newWeight;
}

/********************** Graph  ****************************/

class Vertex;
class Edge;

inline Graph::~Graph() {
    for (Vertex* v : vertexSet) {
        delete v;
    }
}

inline Vertex* Graph::findVertex(long id) const {
    Vertex * v = new Vertex(id);
    auto it = vertexSet.find(v);

    if(it != vertexSet.end()){
        return *it;
    }

    return nullptr;
}

inline bool Graph::addVertex(const Vertex* vertex) {
    if (findVertex(vertex->getId()) != nullptr) {
        return false;  // Vertex with given id already exists
    }
    vertexSet.insert(new Vertex(vertex->getId() , vertex->getLongitude() , vertex->getLatitude()));
    return true;
}

inline bool Graph::removeVertex(long id) {
    auto it = vertexSet.begin();
    while (it != vertexSet.end()) {
        if ((*it)->getId() == id) {
            Vertex* v = *it;
            vertexSet.erase(it);
            delete v;
            return true;
        }
        ++it;
    }
    return false;  // Vertex with given id not found
}

inline bool Graph::addEdge(long sourceId, long destId, double weight) {
    Vertex* source = findVertex(sourceId);
    Vertex* dest = findVertex(destId);
    if (source == nullptr || dest == nullptr) {
        return false;  // Source or destination vertex not found
    }
    source->addAdjEdge(new Edge(source, dest, weight));
    dest->addIncEdge(new Edge(source, dest, weight));
    return true;
}

inline bool Graph::removeEdge(long sourceId, long destId) {
    Vertex* source = findVertex(sourceId);
    Vertex* dest = findVertex(destId);
    if (source == nullptr || dest == nullptr) {
        return false;  // Source or destination vertex not found
    }
    for (Edge* e : source->getAdj()) {
        if (e->getDestination() == dest) {
            source->deleteAdjEdge(e);
            delete e;
            return true;
        }
    }
    return false;  // Edge not found
}

inline bool Graph::addBidirectionalEdge(long sourceId, long destId, double weight) {
    return addEdge(sourceId, destId, weight) && addEdge(destId, sourceId, weight);
}

inline int Graph::getNumVertex() const {
    return vertexSet.size();
}

inline std::unordered_set<Vertex* , VertexHash, VertexEqual> Graph::getVertexSet() const {
    return vertexSet;
}

inline void Graph::clear() {
    this->vertexSet.clear();
    this->edgeSet.clear();
}

inline void Graph::printNodesContente() const {
    for (const auto& vertex : vertexSet) {
        std::cout << "Vertex ID: " << vertex->getId() << std::endl;
        std::cout << "Longitude: " << vertex->getLongitude() << std::endl;
        std::cout << "Latitude: " << vertex->getLatitude() << std::endl;
        std::cout << "Adjacent Edges:" << std::endl;
        for (const auto& edge : vertex->getAdj()) {
            std::cout << "    Destination ID: " << edge->getDestination()->getId() << std::endl;
            std::cout << "    Weight: " << edge->getWeight() << std::endl;
        }
        std::cout << std::endl;
    }
}

/****************** DFS ********************/

/*
void Graph::printMetrics(Metrics metrics) const{
    std::cout << "\033[0;32m" << "-------------------" << "\033[0m";
    std::cout << "\n";
    std::cout << "Average: " << metrics.avg;
    std::cout << "\n";
    std::cout << "Variance: " << metrics.variance;
    std::cout << "\n";
    std::cout << "Maximum Difference: " << metrics.maxDiff;
    std::cout << "\n";
    std::cout << "\033[0;32m" << "-------------------" << "\033[0m";
}*/

#endif /* DA_TP_CLASSES_GRAPH */