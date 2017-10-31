//
//  hGraph.cpp
//  
//
//  Created by Patrick on 10/14/17.
//
//

#include "hGraph.h"


hGraph::hGraph(int size):NUM_NODES(size) {
    _adjMatrix = MatrixXi::Zero(size, size);
    _degVector = Eigen::VectorXi::Zero(size);
    
}

hGraph::hGraph(int size, MatrixXi adjMatrix):NUM_NODES(size) {
    _adjMatrix = MatrixXi::Zero(size, size);
    _adjMatrix = adjMatrix;
    _degVector = Eigen::VectorXi::Zero(size);
    
    for(int i = 0; i < NUM_NODES; i++) {
        for(int j = 0; j < size; j++) {
            _degVector[i] += _adjMatrix(i, j);
        }
    }
}


hGraph::~hGraph() {
    
}

void hGraph::print() {
    
    std::cout << "Adjacency Matrix: " << std::endl << _adjMatrix << std::endl;
    std::cout << "Node Degrees: " << std::endl << _degVector << std::endl;
    
    
}

void hGraph::setMatrix(MatrixXi data) {
}

void hGraph::toFile(std::ofstream &fs) const {
    for(int i = 0; i < NUM_NODES; i++) {
        for(int j = 0; j < NUM_NODES; j++) {
            int temp = _adjMatrix(i,j);
            fs << temp << ",";
            
        }
    }
    
    for(int i = 0; i < NUM_NODES; i++) {
        fs << _degVector[i] << ",";
    }
    
    fs << _eulerChar << "," << _hamiltonian << std::endl;
}

void hGraph::toStream(std::ostream &os) const {

    os << "Adjacency Matrix: " << std::endl << _adjMatrix << std::endl;
    os << "Node Degrees: " << std::endl << _degVector << std::endl;

}

std::ostream &operator << (std::ostream &os, const hGraph &rhs)  {
    rhs.toStream(os);
    return os;
}

std::ofstream &operator << (std::ofstream &fs, const hGraph &rhs)  {
    rhs.toFile(fs);
    return fs;
}









