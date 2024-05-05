//
// Created by jose-costa on 5/5/24.
//

#ifndef PROJETO2_EDGE_H
#define PROJETO2_EDGE_H

#include "Vertex.h"
#include "stdafx.h"

class Edge {
public:
    // Constructor with all parameters
    Edge(const Vertex* source, const Vertex* destination, double weight)
            : source(source), destination(destination), weight(weight) {}


// Setters
    void setSource(const Vertex* newSource) { source = newSource; }
    void setDestination(const Vertex* newDestination) { destination = newDestination; }
    void setWeight(double newWeight) { weight = newWeight; }

    // Getters
    const Vertex* getSource() const { return source; }
    const Vertex* getDestination() const { return destination; }
    double getWeight() const { return weight; }

private:
    const Vertex* source;
    const Vertex* destination;
    double weight;
};



#endif //PROJETO2_EDGE_H
