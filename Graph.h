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
    long getId() const;
    double getLongitude() const;
    double getLatitude() const;
    std::unordered_set<Edge*> getAdj() const;
    void addAdjEdge(Edge* edge);
    void deleteEdge(Edge* edge);

    void setId(long newId);
    void setLongitude(double newLongitude);
    void setLatitude(double newLatitude);
protected:
    long id;
    double longitude;
    double latitude;
    std::unordered_set<Edge*> adj;
    bool visited = false;
    bool processing = false;
    unsigned int indegree = 0;
    double dist = 0;

};

Vertex::Vertex(long id, double longitude, double latitude)
        : id(id), longitude(longitude), latitude(latitude) {}

long Vertex::getId() const {
    return id;
}

double Vertex::getLongitude() const {
    return longitude;
}

double Vertex::getLatitude() const {
    return latitude;
}

std::unordered_set<Edge*> Vertex::getAdj() const {
    return adj;
}

void Vertex::addAdjEdge(Edge* edge) {
    adj.insert(edge);
}

void Vertex::deleteEdge(Edge* edge) {
    adj.erase(edge);
}

void Vertex::setId(long newId) {
    id = newId;
}

void Vertex::setLongitude(double newLongitude) {
    longitude = newLongitude;
}

void Vertex::setLatitude(double newLatitude) {
    latitude = newLatitude;
}

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
    bool addVertex(long id, double longitude, double latitude);
    bool removeVertex(long id);
    bool addEdge(long sourceId, long destId, double weight);
    bool removeEdge(long sourceId, long destId);
    bool addBidirectionalEdge(long sourceId, long destId, double weight);
    int getNumVertex() const;
    std::unordered_set<Vertex*> getVertexSet() const;

private:
    std::unordered_set<Edge*> edgeSet;
    std::unordered_set<Vertex*> vertexSet;
};

/********************** Edge  ****************************/


Edge::Edge(Vertex* source, Vertex* destination, double weight)
        : source(source), destination(destination), weight(weight) {}

Vertex* Edge::getSource() const {
    return source;
}

Vertex* Edge::getDestination() const {
    return destination;
}

double Edge::getWeight() const {
    return weight;
}

void Edge::setSource(Vertex* newSource) {
    source = newSource;
}

void Edge::setDestination(Vertex* newDestination) {
    destination = newDestination;
}

void Edge::setWeight(double newWeight) {
    weight = newWeight;
}

/********************** Graph  ****************************/

class Vertex;
class Edge;

Graph::~Graph() {
    for (Vertex* v : vertexSet) {
        delete v;
    }
}

Vertex* Graph::findVertex(long id) const {
    for (Vertex* v : vertexSet) {
        if (v->getId() == id) {
            return v;
        }
    }
    return nullptr;
}

bool Graph::addVertex(long id, double longitude, double latitude) {
    if (findVertex(id) != nullptr) {
        return false;  // Vertex with given id already exists
    }
    vertexSet.insert(new Vertex(id, longitude, latitude));
    return true;
}

bool Graph::removeVertex(long id) {
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

bool Graph::addEdge(long sourceId, long destId, double weight) {
    Vertex* source = findVertex(sourceId);
    Vertex* dest = findVertex(destId);
    if (source == nullptr || dest == nullptr) {
        return false;  // Source or destination vertex not found
    }
    source->addAdjEdge(new Edge(source, dest, weight));
    return true;
}

bool Graph::removeEdge(long sourceId, long destId) {
    Vertex* source = findVertex(sourceId);
    Vertex* dest = findVertex(destId);
    if (source == nullptr || dest == nullptr) {
        return false;  // Source or destination vertex not found
    }
    for (Edge* e : source->getAdj()) {
        if (e->getDestination() == dest) {
            source->deleteEdge(e);
            delete e;
            return true;
        }
    }
    return false;  // Edge not found
}

bool Graph::addBidirectionalEdge(long sourceId, long destId, double weight) {
    return addEdge(sourceId, destId, weight) && addEdge(destId, sourceId, weight);
}

int Graph::getNumVertex() const {
    return vertexSet.size();
}

std::unordered_set<Vertex*> Graph::getVertexSet() const {
    return vertexSet;
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