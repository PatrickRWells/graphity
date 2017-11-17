//
//  hGraph.cpp
//  
//
//  Created by Patrick on 10/14/17.
//
//

#include "hGraph.h"
#include "absHamiltonian.h"


hGraph::hGraph(int size):NUM_NODES(size) {
    _adjMatrix = MatrixXi::Zero(size, size);
    _degVector = Eigen::VectorXi::Zero(size);
    
}

hGraph::hGraph() {

    _adjMatrix = MatrixXi::Zero(1, 1);
    _degVector = Eigen::VectorXi::Zero(1);
    _hamiltonian = 0.0;

}

hGraph::hGraph(int size, MatrixXi adjMatrix):NUM_NODES(size) {
    _adjMatrix = MatrixXi::Zero(size, size);
    _adjMatrix = adjMatrix;
    _degVector = Eigen::VectorXi::Zero(size);
    _hamiltonian = 0.0;
    
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

void hGraph::setHamiltonian(double val) {
    _hamiltonian = val;
}

void hGraph::setMatrix(MatrixXi data) {
}

void hGraph::setMatrix(int size, MatrixXi data) {
    NUM_NODES = size;
    _adjMatrix.resize(size, size);
    _adjMatrix = data;
    _degVector.resize(size);
    _degVector = Eigen::VectorXi::Zero(size);
    _hamiltonian = 0.0;
    
    for(int i = 0; i < NUM_NODES; i++) {
        for(int j = 0; j < size; j++) {
            _degVector[i] += _adjMatrix(i, j);
        }
    }
    
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
    os << "Value of most recently used hamiltonian: " << _hamiltonian << std::endl;

}

int hGraph::getDegree(int node) {
    if(node < 0 || node >= NUM_NODES) {
        std::cerr << "Critical error, only acceptable node values are between 0 and " << NUM_NODES - 1 << std::endl;
        exit(4);
    }
    
    return _degVector(node);
    
    
}

int hGraph:: getSize() {
    return NUM_NODES;
}

void hGraph::accept(absHamiltonian &ham) {
    std::cout << *this << std::endl;
    ham.calculate(*this);
    
}
    

std::ostream &operator << (std::ostream &os, const hGraph &rhs)  {
    rhs.toStream(os);
    return os;
}

std::ofstream &operator << (std::ofstream &fs, const hGraph &rhs)  {
    rhs.toFile(fs);
    return fs;
}


hGraph * readGraphFile(int &num) { //reads graphs from a CSV, returns a pointer to an array of hGraphs and sets variable num to the number of graphs read from the file
    
    int size;
    char cSize;
    std::string filename;
    std::ifstream input;
    
    while(true) {
        std::cout << "Input graph data filename: "; //gets input file name
        std::cin >> filename;
        input.open(filename);
        if(input.good()) {                  //makes sure file exists
            break;
        }
        else {
            std::cout << "Invalid file. Check that the filename is spelled correctly and that it is in the main folder." << std::endl;
        }
    }
    
    
    std::string line;
    getline(input, line);
    size = line[0] - '0'; //if a character represents an integer, subtracting character 0 from it converts its binary representation to the actual integer.
    
    int linesize = 2*size*size; //figures out how much of the line must be grabbed (including commas) to get the adjacency matrix
    int data[size*size];        //creates array that will temporarily hold adjacency matrix data;
    hGraph * graphData;
    
    int numRead = 0;
    
    while(true) {
        getline(input,line); //gets the next line of the CSV file
        
        if(input.eof()) {  //checks if end of CSV file has been reached;
            break;
        }
        
        numRead++;
    }
    input.clear();
    input.seekg(0, std::ios::beg);
    getline(input, line);
    graphData = new hGraph[numRead];
    for(int i = 0; i < numRead; i++) {
        getline(input,line); //gets the next line of the CSV file
        
        MatrixXi adjMatrix = MatrixXi::Zero(size, size);
        
        for(int j = 0; j < 2*size*size; j += 2) {
            adjMatrix(j/(2*size), (j/2) % size) = line[j] - '0';   //adds the integer to the adjacency matrix in the appropriate position
        }
        
        graphData[i].setMatrix(size, adjMatrix);
        
    }
    
    std::cout << "Number of graphs read " << numRead << std::endl;
    num = numRead;
    
    return graphData;

}







