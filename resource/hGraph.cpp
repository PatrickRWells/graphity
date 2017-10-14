//
//  hGraph.cpp
//  
//
//  Created by Patrick on 10/14/17.
//
//

#include "hGraph.h"

hGraph::hGraph(int size):NUM_NODES(size) {
    _adjMatrix.resize(NUM_NODES, NUM_NODES);
    _degVector.resize(NUM_NODES);
    
}

hGraph::~hGraph() {
    
}

void hGraph::print() {
    
    std::cout << "Adjacency Matrix: " << std::endl << _adjMatrix << std::endl;
    std::cout << "Node Degrees: " << std::endl << _degVector << std::endl;
    
    
}

