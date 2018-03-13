//
//  hGraph.cpp
//  
//
//  Created by Patrick on 10/14/17.
//
//

#include "hGraph.h"
#include "absHamiltonian.h"
/*
------CONTENTS------
 hGraph Class
    -Constructors
    -Destructor
    -Getters
    -Setters
    -Calculaton functions
    -I/O functions
    -hGraph rseource functions
 
 hList class
    -Constructor
    -Destructor
    -Getters
    -List management functions
    -I/O functions
 
 hNode Class
    -Constructor
    -Destructor
    -Node management functions
    -I/O functions
 
 Other resource functions
    -Graph generators
    -I/O  utility functions
    -Matrix mutation functions
 
 */


/////-------------------------BEGIN HGRAPH CLASS IMPLEMENTATION-------------------------/////

//---------------------------CONSTRUCTORS---------------------------//

hGraph::hGraph(int size):NUM_NODES(size) { //Basic constructor. Creates an graph of given size with no connections
    _adjMatrix = MatrixXi::Zero(size, size);
    _degVector = Eigen::VectorXi::Zero(size);
    _numCliques = std::vector <int> (size, 0);
    
}

hGraph::hGraph() {  //default constructor. Creates placeholder. Should not be used often.

    _adjMatrix = MatrixXi::Zero(1, 1);
    _degVector = Eigen::VectorXi::Zero(1);
    _hamiltonian = 0.0;

}

hGraph::hGraph(int size, MatrixXi adjMatrix):NUM_NODES(size) { //Takes a size and a previously defined adjacency matrix and creates an hGraph object.
    _adjMatrix = MatrixXi::Zero(size, size);
    _adjMatrix = adjMatrix;
    _degVector = Eigen::VectorXi::Zero(size);
    _hamiltonian = 0.0;
    _dimension = 0.0; //calcDimension();
    _numCliques = std::vector<int> (size, 0);
    for(int i = 0; i < NUM_NODES; i++) {
        for(int j = 0; j < size; j++) {
            _degVector[i] += _adjMatrix(i, j);
        }
    }
}

//---------------------------END CONSTRUCTORS---------------------------//


//---------------------------DESTRUCTOR---------------------------//


hGraph::~hGraph() { //Since the hGraph class does not include any pointers, there is no work for this constructor to do.
    
}

//---------------------------END DESTRUCTOR---------------------------//

//---------------------------GETTERS---------------------------//

double hGraph::getDimension() { //simple getter function
    
    return _dimension;
    
}

double hGraph::getHam() {  //simple getter function
    return _hamiltonian;
}

int hGraph::getEulerChar() { //simple getter function
    calcEulerChar();
    return _eulerChar;
}

int hGraph:: getSize() { //simple getter function
    return NUM_NODES;
}

void hGraph::numCliques() { //simply outputs the number of cliques of any given size. The cliques are stored in an vector attribute called "_numCliques"
    for (int i = 0; i < NUM_NODES; i++) {
        std::cout << _numCliques.at(i);
    }
    std::cout << std::endl;
}

int hGraph::getDegree(int node) { //simply gets the degree of a given node.
    if(node < 0 || node >= NUM_NODES) {
        std::cerr << "Critical error, only acceptable node values are between 0 and " << NUM_NODES - 1 << std::endl;
        exit(4);
    }
    
    return _degVector(node);
    
    
}

bool hGraph::isConnected(int row, int column) { //checks if two nodes are connected.
    if(_adjMatrix(row, column) == 1) {
        return true;
    }
    
    return false;
    
}



//---------------------------END GETTERS---------------------------//

//---------------------------SETTERS---------------------------//

void hGraph::setHamiltonian(double val) { //sets the value of the hamoiltonian
    _hamiltonian = val;
}

void hGraph::flipEdge(int nodeA, int nodeB) {    
    if(isConnected(nodeA, nodeB)) {
        _adjMatrix(nodeA, nodeB) = 0;
        _adjMatrix(nodeB, nodeA) = 0;
        _degVector[nodeA]--;
        _degVector[nodeB]--;
        
    }
    else {
        _adjMatrix(nodeA, nodeB) = 1;
        _adjMatrix(nodeB, nodeA) = 1;
        _degVector[nodeA]++;
        _degVector[nodeB]++;
    }

    cliquesFound = false;
    
}

//---------------------------END SETTERS---------------------------//


//---------------------------CALCULATIONS---------------------------//

void hGraph::calcEulerChar() { //Based on the definition by Oliver Knills. Requires counting of all cliques in the graph.
    if(!cliquesFound) {         //checks to make sure the cliques have been counted, does so if they have not.
        countCliques();
    }
    int sum = 0;
    for(int i = 0; i < _numCliques.size(); i++) { //See Knill, "On the Dimensionality and Euler Characteristic of Random Graphs"
        if((i+1) % 2 == 1) {                      //for information regarding this algorithm.
            sum -= _numCliques[i];
        }
        else {
            sum += _numCliques[i];
        }
    }
    _eulerChar = sum;
    
}

double hGraph::calcDimension() { //Calculates dimensionality recursively. See Knill, "On the dimensionanity and Euler Characteristic of Random Graphs"
    if(NUM_NODES == 0) {         //for information regarding this algorithm.
        return -1;
    }
    
    else {
        double dimen = 0.0;
        for(int i = 0; i < NUM_NODES; i++ ) {
            dimen += unitSphere(i).getDimension();
        }
        dimen = dimen/(static_cast<double>(NUM_NODES));
        return (dimen + 1);
        
    }
    
    
}

hGraph hGraph::unitSphere(int node) { //outputs the unit sphere around a given. Node the united sphere contains
    //all nodes connected to the given node, as well as all edges between those nodes
    
    
    MatrixXi newAdj(NUM_NODES, NUM_NODES);
    newAdj = _adjMatrix;
    int rem = 0;
    for (int i = 0; i < NUM_NODES; i++) {
        
        if(_adjMatrix(i, node) == 0) {      //Basic algorithm. Simply removes the row and column for a node from the adajecny matrix
            removeColumn(newAdj, i - rem);  //if that node is not connected to the input node.
            removeRow(newAdj, i - rem);
            rem++;
            
            
        }
        
    }
    
    
    hGraph temp(newAdj.rows(), newAdj);
    return temp;
    
}

void hGraph::countCliques() { //Declares necessary objects for initial call of cliqueSearch, then counts all cliques.
    std::vector<int> vectorR(0);
    std::vector<int> vectorP(NUM_NODES);
    for(int i = 0; i < NUM_NODES; i++) {
        vectorP[i] = i;
    }
    
    cliqueSearch(vectorR, vectorP);
    std::cout << "Number of cliques of each size" << std::endl;
    for(int i = 0; i < NUM_NODES; i++) {
        std::cout << i+1 << ": " << _numCliques[i] << std::endl;
        
    }
    cliquesFound = true;
}

void hGraph::cliqueSearch(std::vector<int> R, std::vector<int> P) { //Uses a modified Bron-Kerbosch algorithm to count all ciques of every size.
    
    
    if(R.size() > 0) {
        _numCliques[R.size()-1]++;
    }
    
    if(P.size() > 0) {
        
        while(!P.empty()) {
            std::vector<int> newVectorR(R);
            std::vector<int> newVectorP(0);
            
            int item = P.back();
            newVectorR.push_back(item);
            
            for(int i = P.size() - 1; i >= 0 ; i--) {
                
                if(isConnected(item, P[i])) {
                    newVectorP.push_back(P[i]);
                }
            }
            
            
            cliqueSearch(newVectorR, newVectorP);
            P.pop_back();
            
        }
    }
}


//---------------------------END CALCULATIONS---------------------------//


//---------------------------I/O FUNCTIONS---------------------------//

void hGraph::toFile(std::ofstream &fs) const { //outputs hGraphs to a file in a CSV format. See documentation for information on formatting and usage.
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

std::ofstream &operator << (std::ofstream &fs, const hGraph &rhs)  { //allows use of << operator to output hGraphs to files.
    rhs.toFile(fs);
    return fs;
}


void hGraph::toStream(std::ostream &os) const { //outputs hGraphs in a human-readable format for standard output. See documentation for information on formatting and usage.
    
    os << "Adjacency Matrix: " << std::endl << _adjMatrix << std::endl;
    os << "Node Degrees: " << std::endl << _degVector << std::endl;
    os << "Dimensionality: " << _dimension << std::endl;
    os << "Euler Characteristic: " << _eulerChar << std::endl;
    os << "Value of most recently used hamiltonian: " << _hamiltonian << std::endl;
    
}

std::ostream &operator << (std::ostream &os, const hGraph &rhs)  { //Allows use of << operator to output hGraphs to standard output
    rhs.toStream(os);
    return os;
}




void hGraph::print() { //No longer used.
    
    std::cout << "Adjacency Matrix: " << std::endl << _adjMatrix << std::endl;
    std::cout << "Node Degrees: " << std::endl << _degVector << std::endl;
    
    
}

//---------------------------END I/O FUNCTIONS---------------------------//


//---------------------------HGRAPH RESOURCE FUNCTIONS---------------------------//



void hGraph::setMatrix(int size, MatrixXi data) { //resets the matrix based on new data. Should not be neccesary to use often.
    NUM_NODES = size;
    _adjMatrix.resize(size, size);
    _adjMatrix = data;
    _degVector.resize(size);
    _degVector = Eigen::VectorXi::Zero(size);
    _hamiltonian = 0.0;
    _dimension = calcDimension();
    _numCliques = std::vector<int> (size, 0);
    cliquesFound = false;
    
    
    for(int i = 0; i < NUM_NODES; i++) {
        for(int j = 0; j < size; j++) {
            _degVector[i] += _adjMatrix(i, j);
        }
    }
    
}


void hGraph::accept(absHamiltonian &ham) { //Used in visitor design pattern.
    ham.calculate(*this);
    
}

//---------------------------END HGRAPH RESOURCE FUNCTIONS---------------------------//

/////-------------------------END HGRAPH CLASS IMPLEMENTATION-------------------------/////

    






/////-------------------------BEGIN HLIST CLASS IMPLEMENTATION-------------------------/////

//---------------------------CONSTRUCTOR---------------------------//


hList::hList() { //basic construcotr
    _energy = 0; //all items in the hList should have the same energy.
    _head = nullptr;
}

//---------------------------END CONSTRUCTOR---------------------------//

//---------------------------DESTRUCTOR---------------------------//

hList::~hList() { //clears hList
    if(_head != nullptr) {
        _head->clear();
    }
    delete _head;
}


//---------------------------END DESTRUCTOR---------------------------//


//---------------------------GETTERS---------------------------//

bool hList::isEmpty() {
    return _head == nullptr;
}


double hList::getEnergy() {
    return _energy;
}

//---------------------------END GETTERS---------------------------//

//---------------------------LIST MANAGEMENT FUNCTIONS---------------------------//


void hList::prepend(hGraph * graph) { //adds item to the beginning of the list.
    if(_head == nullptr) {
        _energy = graph->getHam();
    }
    hNode * p = _head;
    hNode * t = new hNode(graph);
    t->_next = p;
    _head = t;
}

void hList::clear() {
    if(_head != nullptr) {
        _head->clear(); //recursively deltetes the items in the list. See implementation of hNode.
        delete _head;
        _head = nullptr;
    }
}

//---------------------------END LIST MANAGEMENT FUNCTIONS---------------------------//

//---------------------------I/O FUNCTIONS---------------------------//



void hList::print() {
    if(_head != nullptr) { //prints items in the list in order. NOTE: should be implemented similarly to file output.
        _head->print();
    }
    else {
        std::cout << "This list is empty" << std::endl;
    }
}

std::ofstream &operator <<(std::ofstream &os, const hList & rhs) { //Allows use of << operator for outputting hLists to files
    rhs.toFile(os);
    return os;
}

void hList::toFile(std::ofstream &fs) const { //outputs items in the list to a file in order.
    if(_head != nullptr) {
        _head->toFile(fs);
    }
}



//---------------------------END I/O FUNCTIONS---------------------------//

/////-------------------------END HLIST CLASS IMPLEMENTATION-------------------------/////


/////-------------------------BEGIN NODE CLASS IMPLEMENTATION-------------------------/////

//---------------------------CONSTRUCTOR---------------------------//


hNode::hNode(hGraph * graph) { //Creates a node with the given graph.
    _data = graph;
    _next = nullptr;
}

//---------------------------END CONSTRUCTOR---------------------------//

//---------------------------DESTRUCTOR---------------------------//



hNode::~hNode() { // Should only be called after clear method.
    _data = nullptr;
    _next = nullptr;

    
}

//---------------------------END DESTRUCTOR---------------------------//

//---------------------------HNODE MANAGEMENT FUNCTIONS---------------------------//

void hNode::clear() { //
    delete _data;
    _data = nullptr;
    if(_next != nullptr) {
        _next->clear();
        delete _next;
        _next = nullptr;
    }
    
}

//---------------------------END HNODE MANAGEMENT FUNCTIONS---------------------------//


//---------------------------I/O FUNCTIONS---------------------------//


void hNode::print() {
    std::cout << *_data << std::endl; //prints out an individual node, then prints the next one if it exists.
    if(_next != nullptr) {
        _next->print();
    }
}

void hNode::toFile(std::ofstream &fs) const { //Outputs hNode data to files. Used within hList class
    fs << *(_data);
    if(_next != nullptr) {
        _next->toFile(fs);
    }
    
}

//---------------------------END I/O FUNCTIONS---------------------------//

/////-------------------------END HNODE CLASS IMPLEMENTATION-------------------------/////



/////-------------------------OTHER RESOURCE FUNCTIONS-------------------------/////

//---------------------------GRAPH GENERATORS---------------------------//

hGraph kGraph(int size) { //generates a complete graph of given size.
    
    MatrixXi adjMatrix = MatrixXi::Zero(size, size);
    for(int i = 0; i < size; i++) {
        for(int j = 0; j < size; j++) {
            if(j != i) {
                adjMatrix(i, j) = 1;
            }
        }
    }
    
    hGraph kGraph(size, adjMatrix);
    return kGraph;
    
    
    
    
}

hGraph randomGraph(int size) { //generates a random graph of a given size.
    
    
    unsigned int max = size*(size-1)/2; //maximum number of edges
    
    std::random_device rd; //c++11 random number generator
    std::mt19937 gen(rd());

    
    
    int fill[max];
    for(int i = 0; i < max; i++ ) { //creates an array that can be mapped to an adjacency matrix and sets all items to zero.
        fill[i] = 0;
    }
    unsigned int fillNum;
    std::uniform_int_distribution<int> dist(max/4, max); //The graph will have at least 1/4th of the maximum number of possible edges
    fillNum = dist(gen);                                 //Ensures graphs are dense enough to be useful.
    std::uniform_int_distribution<int> dist2(0, max-1);  //random number generator that can pick a random edge connection.
    
    
    while(fillNum > 0) { //Randomly fills the array map with 1s representing connections between nodes.
        int i = dist2(gen);
        if (fill[i] == 0){
            fill[i] = 1;
            fillNum--;
            
            
        }
    }
    
    MatrixXi adjMatrix = MatrixXi::Zero(size, size); //creates an adjacency matrix.
    int index = 0;
    for(int k = 1; k < size; k++) {
        for(int m = 0; m < k; m++) {
            adjMatrix(k, m) = fill[index];  //maps the item in the array to its corresponding entry in the bottom half of the adjacency matrix
            adjMatrix(m, k) = fill[index];  //maps the item in the array to its corresponding entry in the top half of the adjacency matrix
            index++;
            
        }
        
    }
    hGraph randGraph(size, adjMatrix);  //Creates an hGraph object using the adjacency matrix.
    return randGraph;
}

//---------------------------END GRAPH GENERATORS---------------------------//


//---------------------------I/O UTILITY FUNCTIONS---------------------------//



hGraph * readGraphFile(int &num) { //reads graphs from a CSV, returns a pointer to an array of hGraphs and sets variable num to the number of graphs read from the file.
                                   //See documentation for how to output graphs to CSV.
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

//---------------------------END I/O UTILITY FUNCTIONS---------------------------//


//---------------------------MATRIX MUTATION FUNCTIONS---------------------------//


void removeRow(Eigen::MatrixXi& matrix, unsigned int rowToRemove) // these two functions are used for generating the unit sphere. Removes rows and columns from a matrix
{
    unsigned int numRows = matrix.rows()-1;
    unsigned int numCols = matrix.cols();
    
    if( rowToRemove < numRows )
        matrix.block(rowToRemove,0,numRows-rowToRemove,numCols) = matrix.block(rowToRemove+1,0,numRows-rowToRemove,numCols);
    
    matrix.conservativeResize(numRows,numCols);
}

void removeColumn(Eigen::MatrixXi& matrix, unsigned int colToRemove)
{
    unsigned int numRows = matrix.rows();
    unsigned int numCols = matrix.cols()-1;
    
    if( colToRemove < numCols )
        matrix.block(0,colToRemove,numRows,numCols-colToRemove) = matrix.block(0,colToRemove+1,numRows,numCols-colToRemove);
    
    matrix.conservativeResize(numRows,numCols);
}

//---------------------------END MATRIX MUTATION FUNCTIONS---------------------------//

/////-------------------------END OTHER RESOURCE FUNCTIONS-------------------------/////


