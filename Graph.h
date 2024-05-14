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
#include <unordered_map>
#include "Logic/Clock.h"

#include "data_structures/MutablePriorityQueue.h"

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

    Vertex(int id, double intitude, double latitude);
    Vertex(int id);
    int getId() const;
    int getIndegree() const;
    int getOutdegree() const;
    double getLongitude() const;
    double getLatitude() const;
    double getDist() const;
    int getNum() const;
    int getLow() const;
    bool isVisited() const;
    bool isProcessing() const;
    Edge * getPath() const;


    std::unordered_map<int, Edge*> getAdj() const;
    std::unordered_map<int, Edge*> getInc() const;
    void addAdjEdge(int edgeId ,Edge* edge);
    void addIncEdge(int edgeId ,Edge* edge);
    void deleteIncEdge(int edgeId);
    void deleteAdjEdge(int edgeId);

    void setDist(double newDist);
    void setId(int newId);
    void setintitude(double newintitude);
    void setLatitude(double newLatitude);
    void setVisited(bool visited);
    void setPath(Edge * e);
    void setOutdegree(int outdegree);
    void setIndegree(int indegree);
    void setNum(int i);
    void setLow(int i);
    void setProcessing(const bool & _processing);

    bool operator<(Vertex & vertex) const;
    bool operator==(Vertex & vertex) const;
    void updateDegrees();
    friend class MutablePriorityQueue<Vertex>;



protected:
    int id;
    double intitude;
    double latitude;
    std::unordered_map<int, Edge*> incomingEdges;
    std::unordered_map<int, Edge*> adjacentEdges;
    bool visited = false;
    bool processing = false;
    unsigned int indegree = 0;
    unsigned int outdegree = 0;
    double dist = 0;
    Edge * path = nullptr;
    int queueIndex = 0;

    int num = 0;
    int low = 0;
};

inline Vertex::Vertex(int id, double intitude, double latitude)
        : id(id), intitude(intitude), latitude(latitude) {}

inline Vertex::Vertex(int id) : id(id) {}

inline int Vertex::getId() const {
    return id;
}

inline double Vertex::getLongitude() const {
    return intitude;
}

inline double Vertex::getLatitude() const {
    return latitude;
}

inline int Vertex::getIndegree() const {
    return indegree;
}

inline int Vertex::getOutdegree() const {
    return outdegree;
}

inline std::unordered_map<int, Edge*> Vertex::getAdj() const {
    return adjacentEdges;
}

inline std::unordered_map<int, Edge*> Vertex::getInc() const {
    return incomingEdges;
}

inline void Vertex::addAdjEdge(int edgeId , Edge* edge) {
    adjacentEdges[edgeId] = edge;
}

inline void Vertex::addIncEdge(int edgeId , Edge* edge) {

    incomingEdges[edgeId] = edge;
}

inline void Vertex::deleteIncEdge(int edgeId) {
    incomingEdges.erase(edgeId);
}

inline void Vertex::deleteAdjEdge(int edgeId) {
    adjacentEdges.erase(edgeId);
}

inline void Vertex::setId(int newId) {
    id = newId;
}

inline void Vertex::setintitude(double newintitude) {
    intitude = newintitude;
}

inline void Vertex::setLatitude(double newLatitude) {
    latitude = newLatitude;
}

inline void Vertex::setIndegree(int indegree) {
   this->indegree = indegree;
}

inline void Vertex::setOutdegree(int outdegree) {
    this->outdegree = outdegree;
}

inline double Vertex::getDist() const {
    return this->dist;
}

inline void Vertex::setDist(double newDist) {
    this->dist = newDist;
}

inline bool Vertex::isVisited() const {
    return this->visited;
}

inline void Vertex::setVisited(bool visited) {
    this->visited = visited;
}

inline bool Vertex::operator<(Vertex &vertex) const {
    return this->dist < vertex.dist;
}

inline Edge *Vertex::getPath() const {
    return this->path;
}

inline void Vertex::setPath(Edge *e) {
    this->path = e;
}

inline bool Vertex::operator==(Vertex &vertex) const {
    return this->id == vertex.getId();
}

inline void Vertex::updateDegrees() {
    indegree = incomingEdges.size();
    outdegree = adjacentEdges.size();
}

inline void Vertex::setNum(int i) {
    this->num = i;
}

inline int Vertex::getNum() const {
    return num;
}

inline int Vertex::getLow() const{
    return low;
}

inline void Vertex::setLow(int i) {
    this->low = i;
}

inline bool Vertex::isProcessing() const {
    return this->processing;
}

inline void Vertex::setProcessing(const bool & _processing) {
    this->processing = _processing;
}

/********************** Edge  ****************************/


class Edge {
public:
    Edge(Vertex* source, Vertex* destination, int edgeId , double weight);
    [[nodiscard]] Vertex* getSource() const;
    [[nodiscard]] Vertex* getDestination() const;
    [[nodiscard]] double getWeight() const;
    [[nodiscard]] double getPheromones() const;
    [[nodiscard]] int getId() const;
    bool isSelected() const;

    void setSource(Vertex* newSource);
    void setDestination(Vertex* newDestination);
    void setWeight(double newWeight);
    void setId(int newId);
    void setSelected(bool selected);
    void setPheromones(double pheromones);
private:
    Vertex* source;
    Vertex* destination;
    double weight;
    int id;
    bool selected;
    double pheromones;
};

/********************** Graph  ****************************/


class Graph {
public:
    ~Graph();
    Edge* findEdge(int sourceId, int destId) const;
    Vertex* findVertex(int id) const;
    bool addVertex(int id , double longitude , double latitude);
    bool removeVertex(int id);
    bool addEdge(int sourceId, int destId,int edgeId, double weight) ;
    bool removeEdge(int sourceId, int destId) const;
    bool addBidirectionalEdge(int sourceId, int destId,int edgeId, double weight) ;
    int getNumVertex() const;
    Clock getClock();
    std::unordered_map<int, Vertex*> getVertexSet() const;
    std::unordered_map<int, Edge*> getEdgeSet() const;

    double Harverstein(double longitude1, double latitude1, double longitude2, double latitude2) const;
    void clear();
    void printNodesContente() const;
    void printGraphInfo() const;
    void populate_in_and_out_degree();
    bool makeFullyConnected() ;
private:

    Clock clock;

    std::unordered_map<int, Edge*> edgeSet;
    std::unordered_map<int, Vertex*> vertexSet;
};

/********************** Edge  ****************************/


inline Edge::Edge(Vertex* source, Vertex* destination, int edgeId , double weight)
        : source(source), destination(destination), id(edgeId) , weight(weight) {}

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

inline int Edge::getId() const {
    return id;
}

inline void Edge::setId(int newId) {
    id = newId;
}

inline bool Edge::isSelected() const {
    return selected;
}

inline void Edge::setSelected(bool newSelected){
    this->selected = newSelected;
}

inline double Edge::getPheromones() const {
    return pheromones;
}

inline void Edge::setPheromones(double pheromones) {
    this->pheromones = pheromones;
}

/********************** Graph  ****************************/

inline Graph::~Graph() {
    for (auto& pair : vertexSet) {
        delete pair.second;
    }
    for (auto& pair : edgeSet) {
        delete pair.second;
    }
}

inline Vertex* Graph::findVertex(int id) const {
    auto it = vertexSet.find(id);

    if(it != vertexSet.end()){
        return it->second;
    }

    return nullptr;
}

inline Edge* Graph::findEdge(int sourceId, int destId) const {
    for (const auto& pair : edgeSet) {
        Edge* edge = pair.second;
        if (edge->getSource()->getId() == sourceId && edge->getDestination()->getId() == destId) {
            return edge;
        }
    }
    return nullptr;  // Edge not found
}

inline bool Graph::addVertex(int id , double longitude , double latitude) {
    Vertex * vertex = new Vertex(id , longitude , latitude);
    if (findVertex(vertex->getId()) != nullptr) {
        return false;  // Vertex with given id already exists
    }
    vertexSet[vertex->getId()] = vertex;
    return true;
}

inline bool Graph::removeVertex(int id) {
    auto it = vertexSet.find(id);
    if (it != vertexSet.end()) {
        delete it->second;
        vertexSet.erase(it);
        return true;
    }
    return false;  // Vertex with given id not found
}

inline bool Graph::addEdge(int sourceId, int destId, int edgeId , double weight)  {
    Vertex* source = findVertex(sourceId);
    Vertex* dest = findVertex(destId);
    if (source == nullptr || dest == nullptr) {
        return false;  // Source or destination vertex not found
    }
    Edge* edge = new Edge(source, dest, edgeId, weight);

    source->addAdjEdge(edgeId ,  edge);
    dest->addIncEdge(edgeId , edge);

    edgeSet[edge->getId()] = edge;
    return true;
}

inline bool Graph::removeEdge(int sourceId, int destId) const {
    Vertex* source = findVertex(sourceId);
    Vertex* dest = findVertex(destId);

    if (source == nullptr || dest == nullptr) {
        return false;  // Source or destination vertex not found
    }

    for (auto it = source->getAdj().begin(); it != source->getAdj().end(); ++it) {
        auto edges = *it;
        if (it->second->getDestination()->getId() == destId) {
            dest->deleteIncEdge(it->second->getId());
            delete it->second;
            source->getAdj().erase(it);
            return true;
        }
    }
    return false;  // Edge not found
}

inline bool Graph::addBidirectionalEdge(int sourceId, int destId, int edgeId , double weight)  {

    Vertex* source = findVertex(sourceId);
    Vertex* dest = findVertex(destId);
    if (source == nullptr || dest == nullptr) {
        return false;  // Source or destination vertex not found
    }

    Edge* edge = new Edge(source, dest, edgeId, weight);
    edgeId++;
    Edge* reveresEdge = new Edge(dest , source , edgeId , weight);

    source->addAdjEdge(edgeId ,  edge);
    dest->addIncEdge(edgeId , edge);

    source->addIncEdge(edgeId , reveresEdge);
    dest->addAdjEdge(edgeId , reveresEdge);

    edgeSet[edge->getId()] = edge;
    edgeSet[reveresEdge->getId()] = reveresEdge;

    return true;
}

inline int Graph::getNumVertex() const {
    return vertexSet.size();
}

inline std::unordered_map<int, Vertex*> Graph::getVertexSet() const {
    return vertexSet;
}



inline void Graph::clear() {
    for (auto& pair : vertexSet) {
        delete pair.second;
    }
    vertexSet.clear();
    for (auto& pair : edgeSet) {
        delete pair.second;
    }
    edgeSet.clear();
}

inline void Graph::printGraphInfo() const {
    std::cout << "Number of nodes: " << vertexSet.size() << std::endl;
    std::cout << "Number of Edges: " << edgeSet.size() << std::endl;
}

inline void Graph::printNodesContente() const {
    for (const auto& pair : vertexSet) {
        Vertex* vertex = pair.second;
        std::cout << "Vertex ID: " << vertex->getId() << std::endl;
        std::cout << "intitude: " << vertex->getLongitude() << std::endl;
        std::cout << "Latitude: " << vertex->getLatitude() << std::endl;
        std::cout << "Adjacent Edges:" << std::endl;
        for (const auto& adjPair : vertex->getAdj()) {
            std::cout << "    Destination ID: " << adjPair.second->getDestination()->getId() << std::endl;
            std::cout << "    Weight: " << adjPair.second->getWeight() << std::endl;
        }
        std::cout << std::endl;
    }
}

inline Clock Graph::getClock() {
    return this->clock;
}

inline bool Graph::makeFullyConnected() {
    std::unordered_map<int, Vertex*> vertices = getVertexSet();

    for(auto v : vertices){
        if(v.second->getLatitude() == DBL_MAX || v.second->getLongitude() == DBL_MAX && v.second->getInc().size() != (getNumVertex()-1) && v.second->getAdj().size() != (getNumVertex()-1)){
            return false;
        }
    }

    // Iterate through all pairs of vertices
    for (auto it1 = vertices.begin(); it1 != vertices.end(); ++it1) {
        for (auto it2 = std::next(it1); it2 != vertices.end(); ++it2) {
            int sourceId = it1->first;
            int destId = it2->first;

            Vertex * v1 = findVertex(sourceId);
            Vertex * v2 = findVertex(destId);
            Edge * e = nullptr;

            for(auto pair : v1->getAdj()){
                Edge * tempE = pair.second;
                if(tempE->getDestination() == v2){
                    e = tempE;
                    continue;
                }
            }

            if(e == nullptr){
                addBidirectionalEdge(sourceId, destId, edgeSet.size() + 1, Harverstein(v1->getLongitude() , v1->getLatitude() , v2->getLongitude() , v2->getLatitude()));
            }
        }
    }
    return true;
}

inline double Graph::Harverstein(double longitude1, double latitude1, double longitude2, double latitude2) const{
    double lon1_rad = longitude1 * M_PI / 180.0;
    double lat1_rad = latitude1 * M_PI / 180.0;
    double lon2_rad = longitude2 * M_PI / 180.0;
    double lat2_rad = latitude2 * M_PI / 180.0;

    double dlon = lon2_rad - lon1_rad;
    double dlat = lat2_rad - lat1_rad;
    double a = pow(sin(dlat / 2), 2) + cos(lat1_rad) * cos(lat2_rad) * pow(sin(dlon / 2), 2);
    double c = 2 * atan2(sqrt(a), sqrt(1 - a));

    const double earth_radius_km = 6371000.0;

    double distance = earth_radius_km * c;
    return distance;
}

inline std::unordered_map<int, Edge *> Graph::getEdgeSet() const {
    return edgeSet;
}

inline void Graph::populate_in_and_out_degree() {
    for (const auto& pair: vertexSet){
        pair.second->setVisited(false);
        pair.second->setIndegree(0);
        pair.second->setOutdegree(0);
    }

    for (const auto& pair: vertexSet){
        for (auto& e: pair.second->getAdj()){
            if (e.second->isSelected()){
                e.second->getSource()->setOutdegree(e.second->getSource()->getOutdegree() + 1);
                e.second->getDestination()->setIndegree(e.second->getDestination()->getIndegree() + 1);
            }
        }
    }
}


#endif /* DA_TP_CLASSES_GRAPH */
