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

double hGraph::getHam() {
    return _hamiltonian;
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

double hGraph::getEulerChar() {
	return _eulerChar;
}

int hGraph:: getSize() {
    return NUM_NODES;
}

bool hGraph::isConnected(int row, int column) {
    if(_adjMatrix(row, column) == 1) {
        return true;
    }
    
    return false;
    
}

void hGraph::accept(absHamiltonian &ham) {
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

hList::hList() {
    _energy = 0;
    _head = nullptr;
}
bool hList::isEmpty() {
    return _head == nullptr;
}

hList::~hList() {
    if(_head != nullptr) {
        _head->clear();
    }
    delete _head;
}

void hList::clear() {
    if(_head != nullptr) {
        _head->clear();
        delete _head;
        _head = nullptr;
    }
}

double hList::getEnergy() {
    return _energy;
}

void hList::prepend(hGraph * graph) {
    if(_head == nullptr) {
        _energy = graph->getHam();
    }
    hNode * p = _head;
    hNode * t = new hNode(graph);
    t->_next = p;
    _head = t;
}

void hList::print() {
    if(_head != nullptr) {
        _head->print();
    }
    else {
        std::cout << "This list is empty" << std::endl;
    }
}

void hList::toFile(std::ofstream &fs) const {
    if(_head != nullptr) {
        _head->toFile(fs);
    }
}
std::ofstream &operator <<(std::ofstream &os, const hList & rhs) {
    rhs.toFile(os);
    return os;
}


hNode::hNode(hGraph * graph) {
    _data = graph;
    _next = nullptr;
}

hNode::~hNode() {
    _data = nullptr;
    _next = nullptr;

    
}

void hNode::clear() {
    _data = nullptr;
    if(_next != nullptr) {
        _next->clear();
        delete _next;
        _next = nullptr;
    }
    
}

void hNode::print() {
    std::cout << *_data << std::endl;
    if(_next != nullptr) {
        _next->print();
    }
}

void hNode::toFile(std::ofstream &fs) const {
    _data->toFile(fs);
    if(_next != nullptr) {
        _next->toFile(fs);
    }
    
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
    size = stoi(line);
    hGraph * graphData;
    int numRead = 0;

    
    int linesize = 2*size*size; //figures out how much of the line must be grabbed (including commas) to get the adjacency matrix
    int data[size*size];        //creates array that will temporarily hold adjacency matrix data;
    
    
    
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
    std::cout << "Reading in graph data..." << std::endl;
    for(int i = 0; i < numRead; i++) {
        getline(input,line); //gets the next line of the CSV file
        
        MatrixXi adjMatrix = MatrixXi::Zero(size, size);
        
        for(int j = 0; j < 2*size*size; j += 2) {
            adjMatrix(j/(2*size), (j/2) % size) = line[j] - '0';   //adds the integer to the adjacency matrix in the appropriate position
        }

        graphData[i].setMatrix(size, adjMatrix);
                
    }
    
    
    input.close();
    std::cout << "Number of graphs read " << numRead << std::endl;
    num = numRead;
    
    return graphData;

}







