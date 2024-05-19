// Original code by Gonçalo Leão
// Updated by LEIC5 Grupo X In 2024-2025

#ifndef DA_TP_CLASSES_GRAPH
#define DA_TP_CLASSES_GRAPH

/**
 * @file Graph.h
 * @brief This file contains the declarations of the classes Vertex, Edge, and Graph used during the project.
 */

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

class Edge;

#define INF std::numeric_limits<double>::max()

/************************* Vertex  **************************/

/**
 * @class Vertex
 * @brief Represents a vertex in the graph.
 */
class Vertex {
public:
    /**
     * @brief Constructs a Vertex with the given ID, longitude, and latitude.
     * @param id Vertex ID.
     * @param intitude Longitude of the vertex.
     * @param latitude Latitude of the vertex.
     */
    Vertex(int id, double intitude, double latitude);

    /**
     * @brief Constructs a Vertex with the given ID.
     * @param id Vertex ID.
     */
    Vertex(int id);

    /**
     * @brief Returns the ID of the vertex.
     * @param none.
     * @return Vertex ID.
     * @complexity O(1)
     */
    int getId() const;

    /**
     * @brief Returns the longitude of the vertex.
     * @param none.
     * @return Longitude.
     * @complexity O(1)
     */
    int getIndegree() const;

    /**
     * @brief Returns the latitude of the vertex.
     * @param none.
     * @return Latitude.
     * @complexity O(1)
     */
    int getOutdegree() const;

    /**
     * @brief Returns the degree of the vertex.
     * @param none.
     * @return Degree.
     * @complexity O(1)
     */
    int getDegree() const;

    /**
     * @brief Returns the longitude of the vertex.
     * @param none.
     * @return Longitude.
     * @complexity O(1)
     */
    double getLongitude() const;

    /**
     * @brief Returns the latitude of the vertex.
     * @param none.
     * @return Latitude.
     * @complexity O(1)
     */
    double getLatitude() const;

    /**
     * @brief Returns the distance of the vertex.
     * @param none.
     * @return Distance.
     * @complexity O(1)
     */
    double getDist() const;

    /**
     * @brief Returns the number of the vertex.
     * @param none.
     * @return Number.
     * @complexity O(1)
     */
    int getNum() const;

    /**
     * @brief Returns the low of the vertex.
     * @param none.
     * @return Low.
     * @complexity O(1)
     */
    int getLow() const;

    /**
     * @brief Returns the tentativas of the vertex.
     * @param none.
     * @return Tentativas.
     * @complexity O(1)
     */
    int getTentativas() const;

    /**
     * @brief Returns the visited status of the vertex.
     * @param none.
     * @return True if the vertex is visited, false otherwise.
     * @complexity O(1)
     */
    bool isVisited() const;

    /**
     * @brief Returns the processing status of the vertex.
     * @param none.
     * @return True if the vertex is processing, false otherwise.
     * @complexity O(1)
     */
    bool isProcessing() const;

    /**
     * @brief Returns the path of the vertex.
     * @param none.
     * @return Path.
     * @complexity O(1)
     */
    Edge * getPath() const;

    /**
     * @brief Returns the adjacent edges of the vertex.
     * @param none.
     * @return Adjacent edges.
     * @complexity O(1)
     */
    std::vector<Edge*> getAdj() const;

    /**
     * @brief Returns the incoming edges of the vertex.
     * @param none.
     * @return Incoming edges.
     * @complexity O(1)
     */
    std::vector<Edge*> getInc() const;

    /**
     * @brief Adds an adjacent edge to the vertex.
     * @param edgeId Edge ID.
     * @param edge Pointer to the edge.
     * @complexity O(1)
     */
    void addAdjEdge(int edgeId ,Edge* edge);

    /**
     * @brief Adds an incoming edge to the vertex.
     * @param edgeId Edge ID.
     * @param edge Pointer to the edge.
     * @complexity O(1)
     */
    void addIncEdge(int edgeId ,Edge* edge);

    /**
     * @brief Deletes an incoming edge from the vertex.
     * @param edgeId Edge ID.
     * @complexity O(E), where E is the number of edges.
     */
    void deleteIncEdge(int edgeId);

    /**
     * @brief Deletes an adjacent edge from the vertex.
     * @param edgeId Edge ID.
     * @complexity O(E), where E is the number of edges.
     */
    void deleteAdjEdge(int edgeId);

    /**
     * @brief Sets the distance of the vertex.
     * @param newDist New distance.
     * @complexity O(1)
     */
    void setDist(double newDist);

    /**
     * @brief Sets the visited status of the vertex.
     * @param visited New visited status.
     * @complexity O(1)
     */
    void setId(int newId);

    /**
     * @brief Sets the longitude of the vertex.
     * @param newLongitude New longitude.
     * @complexity O(1)
     */
    void setintitude(double newintitude);

    /**
     * @brief Sets the latitude of the vertex.
     * @param newLatitude New latitude.
     * @complexity O(1)
     */
    void setLatitude(double newLatitude);

    /**
     * @brief Sets the tentativas of the vertex.
     * @param newTentativas New tentativas.
     * @complexity O(1)
     */
    void setTentativas(int newTentativas);

    /**
     * @brief Sets the adjacent edges of the vertex.
     * @param newAdj New adjacent edges.
     * @complexity O(1).
     */
    void setAdj(std::vector<Edge*> newAdj);

    /**
     * @brief Sets the visited status of the vertex.
     * @param visited New visited status.
     * @complexity O(1)
     */
    void setVisited(bool visited);

    /**
     * @brief Sets the path of the vertex.
     * @param e Pointer to the edge.
     * @complexity O(1)
     */
    void setPath(Edge * e);

    /**
     * @brief Sets the indegree of the vertex.
     * @param indegree New indegree.
     * @complexity O(1)
     */
    void setOutdegree(int outdegree);

    /**
     * @brief Sets the outdegree of the vertex.
     * @param outdegree New outdegree.
     * @complexity O(1)
     */
    void setIndegree(int indegree);

    /**
     * @brief Sets the number of the vertex.
     * @param i New number.
     * @complexity O(1)
     */
    void setNum(int i);

    /**
     * @brief Sets the low of the vertex.
     * @param i New low.
     * @complexity O(1)
     */
    void setLow(int i);

    /**
     * @brief Sets the processing status of the vertex.
     * @param _processing New processing status.
     * @complexity O(1)
     */
    void setProcessing(const bool & _processing);

    /**
     * @brief Updates the indegree and outdegree of the vertex.
     * @param none.
     * @complexity O(1)
     */
    bool operator<(Vertex & vertex) const;

    /**
     * @brief Compares two vertices.
     * @param vertex Vertex to compare.
     * @return True if the vertices are equal, false otherwise.
     * @complexity O(1)
     */
    bool operator==(Vertex & vertex) const;

    /**
     * @brief Orders the adjacent edges of the vertex.
     * @param none.
     * @complexity O(E log E), where E is the number of edges.
     */
    std::vector<Edge*> orderEdges();

    /**
     * @brief Orders the adjacent edges of the vertex.
     * @param none.
     * @complexity O(E log E), where E is the number of edges.
     */
    void updateDegrees();

    /**
     * @brief Orders the adjacent edges of a vertex.
     * @param v Pointer to the vertex.
     * @complexity O(E log E), where E is the number of edges.
     */
    friend class MutablePriorityQueue<Vertex>;

protected:
    int id;                            ///< Vertex ID.
    double intitude;                   ///< Longitude.
    double latitude;                   ///< Latitude.
    std::vector<Edge*> incomingEdges;  ///< Incoming edges.
    std::vector<Edge*> adjacentEdges;  ///< Adjacent edges.
    bool visited = false;              ///< Visit status.
    bool processing = false;           ///< Processing status.
    unsigned int indegree = 0;         ///< Indegree.
    unsigned int outdegree = 0;        ///< Outdegree.
    double dist = 0;                   ///< Distance.
    Edge* path = nullptr;              ///< Path edge.
    int queueIndex = 0;                ///< Queue index.

    int tentativas = 0;                ///< Tentativas.
    int num = 0;                       ///< Num.
    int low = 0;                       ///< Low.
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

inline std::vector<Edge*> Vertex::getAdj() const {
    return adjacentEdges;
}

inline std::vector<Edge*> Vertex::getInc() const {
    return incomingEdges;
}

inline void Vertex::addAdjEdge(int edgeId , Edge* edge) {
    adjacentEdges.push_back(edge);
}

inline void Vertex::addIncEdge(int edgeId , Edge* edge) {
    incomingEdges.push_back(edge);
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

inline int Vertex::getDegree() const {
    return indegree + outdegree;
}

inline void Vertex::setAdj(std::vector<Edge*> newAdj) {
    this->adjacentEdges = std::move(newAdj);
}

inline int Vertex::getTentativas() const {
    return tentativas;
}

inline void Vertex::setTentativas(int newTentativas) {
    this->tentativas = newTentativas;
}

/********************** Edge  ****************************/

/**
 * @class Edge
 * @brief Represents an edge in the graph.
 */
class Edge {
public:
    /**
     * @brief Constructs an Edge with the given source, destination, ID, and weight.
     * @param source Source vertex.
     * @param destination Destination vertex.
     * @param edgeId Edge ID.
     * @param weight Weight of the edge.
     */
    Edge(Vertex* source, Vertex* destination, int edgeId , double weight);

    /**
     * @brief Returns the source vertex.
     * @param none.
     * @return Pointer to the source vertex.
     * @complexity O(1)
     */
    [[nodiscard]] Vertex* getSource() const;

    /**
     * @brief Returns the destination vertex.
     * @param none.
     * @return Pointer to the destination vertex.
     * @complexity O(1)
     */
    [[nodiscard]] Vertex* getDestination() const;

    /**
     * @brief Returns the weight of the edge.
     * @param none.
     * @return Weight of the edge.
     * @complexity O(1)
     */
    [[nodiscard]] double getWeight() const;

    /**
     * @brief Returns the pheromones of the edge.
     * @param none.
     * @return Pheromones of the edge.
     * @complexity O(1)
     */
    [[nodiscard]] double getPheromones() const;

    /**
     * @brief Returns the ID of the edge.
     * @param none.
     * @return ID of the edge.
     * @complexity O(1)
     */
    [[nodiscard]] int getId() const;

    /**
     * @brief Returns the selection status of the edge.
     * @param none.
     * @return True if the edge is selected, false otherwise.
     * @complexity O(1)
     */
    bool isSelected() const;

    /**
     * @brief Sets the source vertex.
     * @param newSource Pointer to the new source vertex.
     * @complexity O(1)
     */
    void setSource(Vertex* newSource);

    /**
     * @brief Sets the destination vertex.
     * @param newDestination Pointer to the new destination vertex.
     * @complexity O(1)
     */
    void setDestination(Vertex* newDestination);

    /**
     * @brief Sets the weight of the edge.
     * @param newWeight New weight of the edge.
     * @complexity O(1)
     */
    void setWeight(double newWeight);

    /**
     * @brief Sets the ID of the edge.
     * @param newId New ID of the edge.
     * @complexity O(1)
     */
    void setId(int newId);

    /**
     * @brief Sets the selection status of the edge.
     * @param selected New selection status.
     * @complexity O(1)
     */
    void setSelected(bool selected);

    /**
     * @brief Sets the pheromones of the edge.
     * @param pheromones New pheromones of the edge.
     * @complexity O(1)
     */
    void setPheromones(double pheromones);
private:
    Vertex* source;         ///< Source vertex.
    Vertex* destination;    ///< Destination vertex.
    double weight;          ///< Weight.
    int id;                 ///< Edge ID.
    bool selected;          ///< Selection status.
    double pheromones;      ///< Pheromones.
};

/********************** Graph  ****************************/

/**
 * @class Graph
 * @brief Represents a graph consisting of vertices and edges.
 */
class Graph {
public:
    /**
     * @brief Destructor to clean up the graph.
     */
    ~Graph();

    /**
     * @brief Finds a vertex with the given ID.
     * @param id Vertex ID.
     * @complexity O(E) , where E is the number of edges.
     * @return Pointer to the vertex if found, nullptr otherwise.
     */
    Edge* findEdge(int sourceId, int destId) const;

    /**
     * @brief Finds a vertex with the given ID.
     * @param id Vertex ID.
     * @complexity O(1)
     * @return Pointer to the vertex if found, nullptr otherwise.
     */
    Vertex* findVertex(int id) const;

    /**
     * @brief Adds a vertex to the graph.
     * @param id Vertex ID.
     * @param longitude Longitude of the vertex.
     * @param latitude Latitude of the vertex.
     * @complexity O(1)
     * @return True if the vertex was added, false otherwise.
     */
    bool addVertex(int id , double longitude , double latitude);

    /**
     * @brief Removes a vertex from the graph.
     * @param id Vertex ID.
     * @complexity O(E) , where E is the number of edges.
     * @return True if the vertex was removed, false otherwise.
     */
    bool removeVertex(int id);

    /**
     * @brief Adds an edge to the graph.
     * @param sourceId Source vertex ID.
     * @param destId Destination vertex ID.
     * @param edgeId Edge ID.
     * @param weight Weight of the edge.
     * @complexity O(1)
     * @return True if the edge was added, false otherwise.
     */
    bool addEdge(int sourceId, int destId,int edgeId, double weight) ;

    /**
     * @brief Removes an edge from the graph.
     * @param sourceId Source vertex ID.
     * @param destId Destination vertex ID.
     * @complexity O(E) , where E is the number of edges.
     * @return True if the edge was removed, false otherwise.
     */
    bool removeEdge(int sourceId, int destId) const;

    /**
     * @brief Adds a bidirectional edge to the graph.
     * @param sourceId Source vertex ID.
     * @param destId Destination vertex ID.
     * @param edgeId Edge ID.
     * @param weight Weight of the edge.
     * @complexity O(1)
     * @return True if the edge was added, false otherwise.
     */
    bool addBidirectionalEdge(int sourceId, int destId,int edgeId, double weight) ;

    /**
     * @brief Returns the size of the VertexSet.
     * @param none.
     * @complexity O(1).
     * @return an integer holding the size of the vertexSet.
     */
    int getNumVertex() const;

    /**
     * @brief Returns the clock.
     * @param none.
     * @complexity O(1).
     * @return the clock.
     */
    Clock getClock();

    /**
     * @brief Returns the vertex set.
     * @param none.
     * @complexity O(1).
     * @return the vertex set.
     */
    std::unordered_map<int, Vertex*> getVertexSet() const;

    /**
     * @brief Returns the edge set.
     * @param none.
     * @complexity O(1).
     * @return the edge set.
     */
    std::unordered_map<int, Edge*> getEdgeSet() const;

    /**
     * @brief Returns Harverstein distance between two nodes.
     * @param longitude1 the longitude of the first node.
     * @param latitude1 the latitude of the first node.
     * @param longitude2 the longitude of the second node.
     * @param latitude2 the latitude of the second node.
     * @complexity O(logN) , where N is the number whose square root is calculated.
     * @return the edge set.
     */
    double Harverstein(double longitude1, double latitude1, double longitude2, double latitude2) const;

    /**
     * @brief Clears the graph.
     * @param none.
     * @complexity O(V + E), where V is the number of vertices and E is the number of edges.
     */
    void clear();

    /**
     * @brief Prints the content of the nodes.
     * @param none.
     * @complexity O(V + E), where V is the number of vertices and E is the number of edges.
     */
    void printNodesContente() const;

    /**
     * @brief Prints the graph information.
     * @param none.
     * @complexity O(1).
     */
    void printGraphInfo() const;

    /**
     * @brief Populates the indegree and outdegree of the vertices.
     * @param none.
     * @complexity O(V + E), where V is the number of vertices and E is the number of edges.
     */
    void populate_in_and_out_degree();

    /**
     * @brief Orders the adjacent edges of a vertex.
     * @param v Pointer to the vertex.
     * @complexity O(V²), where V is the number of vertexes.
     */
    bool makeFullyConnected() ;

    /**
     * @brief Orders the adjacent edges of a vertex.
     * @param v Pointer to the vertex.
     * @complexity O(E log E), where E is the number of edges.
     */
    void orderVertexAdj(Vertex * v);

private:
    std::unordered_map<int, Vertex*> vertexSet;   ///< Vertex set.
    std::unordered_map<int, Edge*> edgeSet;       ///< Edge set.
    Clock clock;                                  ///< Clock.
};

/********************** Edge  ****************************/

/**
 *
 * @param source
 * @param destination
 * @param edgeId
 * @param weight
 */
inline Edge::Edge(Vertex* source, Vertex* destination, int edgeId , double weight)
        : source(source), destination(destination), id(edgeId) , weight(weight) {}

/**
 * @brief returns the source.
 * @param none.
 * @return the source of the edge.
 * @complexity O(1).
 */

inline Vertex* Edge::getSource() const {
    return source;
}

/**
 * @brief returns the destination.
 * @param none.
 * @return the destination of the edge.
 * @complexity O(1).
 */
inline Vertex* Edge::getDestination() const {
    return destination;
}

/**
 * @brief returns the weight.
 * @param none.
 * @return the weight of the edge.
 * @complexity O(1).
 */
inline double Edge::getWeight() const {
    return weight;
}

/**
 * @brief sets the source.
 * @param newSource
 * @complexity O(1).
 */
inline void Edge::setSource(Vertex* newSource) {
    source = newSource;
}

/**
 * @brief sets the destination.
 * @param newDestination
 * @complexity O(1).
 */
inline void Edge::setDestination(Vertex* newDestination) {
    destination = newDestination;
}

/**
 * @brief sets the weight.
 * @param newWeight
 * @complexity O(1).
 */
inline void Edge::setWeight(double newWeight) {
    weight = newWeight;
}

/**
 * @brief sets the weight.
 * @param newWeight
 * @complexity O(1).
 */

inline int Edge::getId() const {
    return id;
}

/**
 * @brief sets the weight.
 * @param newWeight
 * @complexity O(1).
 */
inline void Edge::setId(int newId) {
    id = newId;
}

/**
 * @brief returns the selection status.
 * @param none.
 * @return the selection status of the edge.
 * @complexity O(1).
 */
inline bool Edge::isSelected() const {
    return selected;
}

/**
 * @brief sets the selection status.
 * @param newSelected
 * @complexity O(1).
 */
inline void Edge::setSelected(bool newSelected){
    this->selected = newSelected;
}

/**
 * @brief returns the pheromones.
 * @param none.
 * @return the pheromones of the edge.
 * @complexity O(1).
 */
inline double Edge::getPheromones() const {
    return pheromones;
}

/**
 * @brief sets the pheromones.
 * @param pheromones
 * @complexity O(1).
 */
inline void Edge::setPheromones(double pheromones) {
    this->pheromones = pheromones;
}

/********************** Graph  ****************************/

/**
 * @brief Destructor to clean up the graph.
 * @param none.
 * @complexity O(V + E), where V is the number of vertices and E is the number of edges.
 */
inline Graph::~Graph() {
    for (auto& pair : vertexSet) {
        delete pair.second;
    }
    for (auto& pair : edgeSet) {
        delete pair.second;
    }
}

/**
 * @brief Finds a vertex with the given ID.
 * @param id Vertex ID.
 * @complexity O(1).
 * @return Pointer to the vertex if found, nullptr otherwise.
 */
inline Vertex* Graph::findVertex(int id) const {
    auto it = vertexSet.find(id);

    if(it != vertexSet.end()){
        return it->second;
    }

    return nullptr;
}

/**
 * @brief Finds an edge between two vertices.
 * @param sourceId Source vertex ID.
 * @param destId Destination vertex ID.
 * @complexity O(E), where E is the number of edges.
 * @return Pointer to the edge if found, nullptr otherwise.
 */
inline Edge* Graph::findEdge(int sourceId, int destId) const {
    for (const auto& pair : edgeSet) {
        Edge* edge = pair.second;
        if (edge->getSource()->getId() == sourceId && edge->getDestination()->getId() == destId) {
            return edge;
        }
    }
    return nullptr;  // Edge not found
}

/**
 * @brief Adds a vertex to the graph.
 * @param id Vertex ID.
 * @param longitude Longitude of the vertex.
 * @param latitude Latitude of the vertex.
 * @complexity O(1).
 * @return True if the vertex was added, false otherwise.
 */
inline bool Graph::addVertex(int id , double longitude , double latitude) {
    Vertex * vertex = new Vertex(id , longitude , latitude);
    if (findVertex(vertex->getId()) != nullptr) {
        return false;  // Vertex with given id already exists
    }
    vertexSet[vertex->getId()] = vertex;
    return true;
}

/**
 * @brief Removes a vertex from the graph.
 * @param id Vertex ID.
 * @complexity O(1).
 * @return True if the vertex was removed, false otherwise.
 */
inline bool Graph::removeVertex(int id) {
    auto it = vertexSet.find(id);
    if (it != vertexSet.end()) {
        delete it->second;
        vertexSet.erase(it);
        return true;
    }
    return false;  // Vertex with given id not found
}

/**
 * @brief Adds an edge to the graph.
 * @param sourceId Source vertex ID.
 * @param destId Destination vertex ID.
 * @param edgeId Edge ID.
 * @param weight Weight of the edge.
 * @complexity O(1).
 * @return True if the edge was added, false otherwise.
 */
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

/**
 * @brief Removes an edge from the graph.
 * @param sourceId Source vertex ID.
 * @param destId Destination vertex ID.
 * @complexity O(E), where E is the number of edges.
 * @return True if the edge was removed, false otherwise.
 */
inline bool Graph::removeEdge(int sourceId, int destId) const {
    Vertex* source = findVertex(sourceId);
    Vertex* dest = findVertex(destId);

    if (source == nullptr || dest == nullptr) {
        return false;  // Source or destination vertex not found
    }

    for (auto it = source->getAdj().begin(); it != source->getAdj().end(); ++it) {
        Edge * edge = *it;

        if (edge->getSource()->getId() == sourceId && edge->getDestination()->getId() == destId) {

            for(auto it2 = dest->getInc().begin() ; it2 != dest->getInc().end() ; it2++){
                Edge * e = *it2;
                if(e->getId() == edge->getId()){
                    it2 = dest->getInc().erase(it2);
                    break;
                }
            }

            for(auto it2 = dest->getAdj().begin() ; it2 != dest->getInc().end() ; it2++){
                Edge * e = *it2;
                if(e->getId() == edge->getId()){
                    it2 = dest->getInc().erase(it2);
                    break;
                }
            }

            for(auto it2 = source->getInc().begin() ; it2 != source->getInc().end() ; it2++){
                Edge * e = *it2;
                if(e->getId() == edge->getId()){
                    it2 = dest->getInc().erase(it2);
                    break;
                }
            }

            for(auto it2 = source->getAdj().begin() ; it2 != source->getInc().end() ; it2++){
                Edge * e = *it2;
                if(e->getId() == edge->getId()){
                    it2 = dest->getInc().erase(it2);
                    break;
                }
            }

            return true;
        }
    }
    return false;  // Edge not found
}

/**
 * @brief Adds a bidirectional edge to the graph.
 * @param sourceId Source vertex ID.
 * @param destId Destination vertex ID.
 * @param edgeId Edge ID.
 * @param weight Weight of the edge.
 * @complexity O(1).
 * @return True if the edge was added, false otherwise.
 */
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

/**
 * @brief Returns the size of the VertexSet.
 * @param none.
 * @complexity O(1).
 * @return an integer holding the size of the vertexSet.
 */
inline int Graph::getNumVertex() const {
    return vertexSet.size();
}

/**
 * @brief Finds a vertex with the given ID.
 * @param id Vertex ID.
 * @complexity O(1).
 * @return Pointer to the vertex if found, nullptr otherwise.
 */
inline std::unordered_map<int, Vertex*> Graph::getVertexSet() const {
    return vertexSet;
}

/**
 * @brief Clears the graph.
 * @param none.
 * @complexity O(V + E), where V is the number of vertices and E is the number of edges.
 */

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

/**
 * @brief Print the Number of nodes and Number of Edges .
 * @param none.
 * @complexity O(1).
 */
inline void Graph::printGraphInfo() const {
    std::cout << "Number of nodes: " << vertexSet.size() << std::endl;
    std::cout << "Number of Edges: " << edgeSet.size() << std::endl;
}

/**
 * @brief Prints the content of the nodes.
 * @param none.
 * @complexity O(V + E), where V is the number of vertices and E is the number of edges.
 */
inline void Graph::printNodesContente() const {
    for (const auto& pair : vertexSet) {
        Vertex* vertex = pair.second;
        std::cout << "Vertex ID: " << vertex->getId() << std::endl;

        if(vertex->getLongitude() == DBL_MAX){
            std::cout << "Longitude: Not defined" << std::endl;
        }
        else{
            std::cout << "Longitude: " << vertex->getLongitude() << std::endl;
        }

        if(vertex->getLatitude() == DBL_MAX){
            std::cout << "Longitude: Not defined" << std::endl;
        }
        else{
            std::cout << "Longitude: " << vertex->getLongitude() << std::endl;
        }

        std::cout << "Adjacent Edges:" << std::endl;
        for (const auto& adjPair : vertex->getAdj()) {
            std::cout << "    Destination ID: " << adjPair->getDestination()->getId() << std::endl;
            std::cout << "    Weight: " << adjPair->getWeight() << std::endl;
        }
        std::cout << std::endl;
    }
}


/**
 * @brief Returns the clock.
 * @param none.
 * @complexity O(1).
 * @return the clock.
 */
inline Clock Graph::getClock() {
    return this->clock;
}

/**
     * @brief Orders the adjacent edges of a vertex.
     * @param v Pointer to the vertex.
     * @complexity O(V²), where V is the number of vertexes.
     */
inline bool Graph::makeFullyConnected() {
    // Get all vertices and edges
    std::unordered_map<int, Vertex*> vertices = getVertexSet();
    std::unordered_map<int, Edge*> edges = getEdgeSet();

    // Check if the graph is already fully connected
    int numVertices = getNumVertex();
    int expectedEdges = numVertices * (numVertices - 1);
    if (edges.size() == expectedEdges || edges.size() == expectedEdges / 2) {
        return false;
    }

    // Check if all vertices have defined coordinates
    for (const auto& v : vertices) {
        if (v.second->getLatitude() == DBL_MAX || v.second->getLongitude() == DBL_MAX) {
            return false;
        }
    }

    // Create a set to store unique pairs of connected vertices
    std::vector<std::vector<bool>> connections(numVertices, std::vector<bool>(numVertices, false));
    for (const auto& e : edges) {
        int sourceId = e.second->getSource()->getId();
        int destId = e.second->getDestination()->getId();
        connections[sourceId][destId] = true;
        connections[destId][sourceId] = true; // Because of bidirectional edges
    }

    // Iterate through all pairs of vertices to add missing edges
    for (auto it1 = vertices.begin(); it1 != vertices.end(); ++it1) {
        for (auto it2 = std::next(it1); it2 != vertices.end(); ++it2) {
            int sourceId = it1->second->getId();
            int destId = it2->second->getId();

            // If there is no edge between sourceId and destId, add it
            if(!connections[sourceId][destId]){
                Vertex* v1 = it1->second;
                Vertex* v2 = it2->second;
                double weight = Harverstein(v1->getLongitude(), v1->getLatitude(), v2->getLongitude(), v2->getLatitude());
                addBidirectionalEdge(sourceId, destId, edgeSet.size() + 1, weight);
                connections[sourceId][destId] = true;
                connections[destId][sourceId] = true;
            }
        }
    }
    return true;
}

/**
     * @brief Returns Harverstein distance between two nodes.
     * @param longitude1 the longitude of the first node.
     * @param latitude1 the latitude of the first node.
     * @param longitude2 the longitude of the second node.
     * @param latitude2 the latitude of the second node.
     * @complexity O(logN) , where N is the number whose square root is calculated.
     * @return the edge set.
     */
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

/**
 * @brief Returns the edge set.
 * @param none.
 * @complexity O(1).
 * @return the edge set.
 */
inline std::unordered_map<int, Edge *> Graph::getEdgeSet() const {
    return edgeSet;
}

/**
 * @brief Populates the indegree and outdegree of the vertices.
 * @param none.
 * @complexity O(V + E), where V is the number of vertices and E is the number of edges.
 */
inline void Graph::populate_in_and_out_degree() {
    for (const auto& pair: vertexSet){
        pair.second->setVisited(false);
        pair.second->setIndegree(0);
        pair.second->setOutdegree(0);
    }

    for (const auto& pair: vertexSet){
        for (Edge * e: pair.second->getAdj()){
            if (e->isSelected()){
                e->getSource()->setOutdegree(e->getSource()->getOutdegree() + 1);
                e->getDestination()->setIndegree(e->getDestination()->getIndegree() + 1);
            }
        }
    }
}

/**
 * @brief Orders the adjacent edges of a vertex.
 * @param v Pointer to the vertex.
 * @complexity O(E log E), where E is the number of edges.
 */
inline void Graph::orderVertexAdj(Vertex * v){
    auto comp = [] (Edge * e1, Edge * e2){
        if (e1->getDestination()->getDegree() < e2->getDestination()->getDegree()) {
            return true;
        }
        else if (e1->getDestination()->getDegree() == e2->getDestination()->getDegree() and e1->getWeight() < e2->getWeight()){
            return true;
        }
        return false;
    };

    std::sort(v->getAdj().begin(), v->getAdj().end(), comp);
}


#endif /* DA_TP_CLASSES_GRAPH */
