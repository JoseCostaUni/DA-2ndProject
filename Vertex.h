//
// Created by jose-costa on 5/5/24.
//

#ifndef PROJETO2_VERTEX_H
#define PROJETO2_VERTEX_H

#include "Edge.h"
#include "stdafx.h"

class Vertex {
public:
    // Constructors
    Vertex(long id, double longitude, double latitude, const std::vector<Edge*>& edges, const std::vector<Edge*>& incomingEdges)
            : id(id), longitude(longitude), latitude(latitude), adjEdges(edges), incomingEdges(incomingEdges) {}

    Vertex(long id) : id(id), longitude(0.0), latitude(0.0) {} // Default values for longitude and latitude

    // Setters
    void setId(long newId) { id = newId; }
    void setLongitude(double newLongitude) { longitude = newLongitude; }
    void setLatitude(double newLatitude) { latitude = newLatitude; }
    void setAdjacentEdges(const std::vector<Edge*>& newEdges) { adjEdges = newEdges; }
    void setIncomingEdges(const std::vector<Edge*>& newEdges) { incomingEdges = newEdges; }

    // Getters
    long getId() const { return id; }
    double getLongitude() const { return longitude; }
    double getLatitude() const { return latitude; }
    const std::vector<Edge*>& getAdjacentEdges() const { return adjEdges; }
    const std::vector<Edge*>& getIncomingEdges() const { return incomingEdges; }

    bool addAdjEdge(Edge* edge) {
        adjEdges.push_back(edge);
        return true; // Assuming always successful for simplicity
    }

    bool addIncomingEdge(Edge* incomingEdge) {
        incomingEdges.push_back(incomingEdge);
        return true; // Assuming always successful for simplicity
    }

    void removeEdge()

    // Hash function for unordered_set
    struct Hash {
        size_t operator()(const Vertex& vertex) const {
            return std::hash<long>()(vertex.id);
        }
    };

private:
    long id;
    double longitude;
    double latitude;
    std::vector<Edge*> adjEdges;
    std::vector<Edge*> incomingEdges;
};

#endif //PROJETO2_VERTEX_H
