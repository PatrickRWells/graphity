    //
//  hGraph.cpp
//  
//
//  Created by Patrick on 10/14/17.
//
//

#include <algorithm>
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
    -Other Calculations
 
 */


/////-------------------------BEGIN HGRAPH CLASS IMPLEMENTATION-------------------------/////

//---------------------------CONSTRUCTORS---------------------------//

hGraph::hGraph(int size):NUM_NODES(size) { //Basic constructor. Creates an graph of given size with no connections
    _adjMatrix = MatrixXi::Zero(size, size);
    _degVector = Eigen::VectorXi::Zero(size);
    _numCliques = std::vector <int> (size, 0);
    for(int i = 0; i < size; i++) {
        _eccentricities.push_back(0);
        
    }
    
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
    _dimension = 0; //this is not initialized by default due to the high computing cost. ???Why is there a phase transition at 36 nodes???
    _numCliques = std::vector<int> (size, 0);
    for(int i = 0; i < NUM_NODES; i++) {
        _eccentricities.push_back(0);
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
    if(_dimension == 0) {
        calcDimension();
    }
    return _dimension;
    
}

std::vector<double> hGraph::getSpectralDimen() {
    if (_spectralDimen.size() == 0) {
        calcSpectralDimen();
    }
    return _spectralDimen;
    
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

std::vector<int> hGraph::getEccentricity() {
    if(_eccentricities[0] == 0 && _eccentricities[NUM_NODES/2] == 0) {
        calculateEccen();
    }
    return _eccentricities;
    
}

int hGraph::getDiameter() {
    std::vector<int> temp = getEccentricity();
    int val = 0;
    for(int i = 0; i < temp.size(); i++ ) {
        if(temp[i] > val) {
            val = temp[i];
        }
    }
    
    return val;
    
}

std::vector<double> hGraph::getHausdorffDimen() {
    
    if(_hausdorffDimen.size() == 0) {
        calculateHausDimen();
        
    }
    
    return _hausdorffDimen;
    
}



//---------------------------END GETTERS---------------------------//

//---------------------------SETTERS---------------------------//

void hGraph::setHamiltonian(double val) { //sets the value of the hamoiltonian
    _hamiltonian = val;
}

void hGraph::setThreads(int threads) {
    _numThreads = threads;
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
    _dimension = 0;
    
}

void hGraph::flipEdge(std::vector<int> nodeA, std::vector<int> nodeB) {
    for(int i = 0; i < nodeA.size(); i++) {
        if(isConnected(nodeA[i], nodeB[i])) {
            _adjMatrix(nodeA[i], nodeB[i]) = 0;
            _adjMatrix(nodeB[i], nodeA[i]) = 0;
            _degVector[nodeA[i]]--;
            _degVector[nodeB[i]]--;
        
        }
        else {
            _adjMatrix(nodeA[i], nodeB[i]) = 1;
            _adjMatrix(nodeB[i], nodeA[i]) = 1;
            _degVector[nodeA[i]]++;
            _degVector[nodeB[i]]++;
        }
    }
    
    
}

void hGraph::acceptPartial(double partial) {
    _hamiltonian += partial;
    
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
            sum += _numCliques[i];
        }
        else {
            sum -= _numCliques[i];
        }
    }
    _eulerChar = sum;
    
}

MatrixXi hGraph::getShortestPaths() {
    MatrixXi pathMatrix(NUM_NODES, NUM_NODES);
    pathMatrix = _adjMatrix;
    for(int i = 0; i < NUM_NODES; i++) {
        for (int j = 0; j < NUM_NODES; j++ ) {
            
            if((i != j) && (pathMatrix(i,j) == 0)) {
                pathMatrix(i,j) = 10000000;
            }
            
        }
    }
    
    
    for(int k = 0; k < NUM_NODES; k++) {
        for(int i = 0; i < NUM_NODES; i++) {
            for(int j = 0; j < NUM_NODES; j++) {
                
                if(pathMatrix(i,j) > (pathMatrix(i,k) + pathMatrix(j, k))) {
                    pathMatrix(i,j) = (pathMatrix(i,k) + pathMatrix(j,k));
                }
                
                
            }
        }
    }
    return pathMatrix;
    
    
}

hGraph hGraph::compliment() {
    hGraph temp(NUM_NODES, _adjMatrix);
    for(int i = 0; i < NUM_NODES; i++) {
        for(int j = i+1; j < NUM_NODES; j++) {
            temp.flipEdge(i, j);
        }
    }

    return temp;
    
    
}

double hGraph::avgDegree() {
    double avg = 0;
    for(int i = 0; i < NUM_NODES; i++) {
        avg += getDegree(i);
    }
    return avg/NUM_NODES;
}


void hGraph::calcDimension() {
    int num = 0;
    int * numAvoided = new int[NUM_NODES];
    int maxEdges = NUM_NODES*(NUM_NODES-1) / 2;
    int edges = 0;
    for(int i = 0; i < NUM_NODES; i++) {
        edges += getDegree(i);
    }
    edges /= 2;
    if(edges == maxEdges) {
        
        std::cout << NUM_NODES - 1 << std::endl;
        return;
    }
    
    for(int i = 0; i < NUM_NODES; i++) {
        numAvoided[i] = 0;
    }
    double *** values;
    values = new double ** [NUM_NODES];
    for(int i = 0; i < NUM_NODES; i++) {
        values[i] = new double * [NUM_NODES];
        for(int j = 0; j < NUM_NODES; j++) {
            values[i][j] = new double[NUM_NODES];
            for(int k = 0; k < NUM_NODES; k++) {
                values[i][j][k] = -2.0;
            }
        }
    }

    double val = 0.0;
    for(int i = 0; i < NUM_NODES; i++) {
        val += dimension(::unitSphere(_adjMatrix, i), values, numAvoided);
        
    }
    val = val/NUM_NODES;
    val++;
    std::cout << val << std::endl;
    //_dimension = val;
    /*for(int i = 0; i < NUM_NODES; i++) {
        
        std::cout << numAvoided[i] << " recursive calls avoided at level " << i << std::endl;

    }*/
}

double hGraph::dimension(MatrixXi amat, double *** data, int  * numAvoided)  {
    double retVal = 0;
    int recDepth = 0;
    int numNodes = 0;
    int numEdges = 0;
    std::vector<int> nodesChecked;
    for(int i = 0; i < NUM_NODES; i++) {
        if(amat(i, 0) == -2) {
            recDepth++;
            nodesChecked.push_back(i);
        }
        else if(amat(i, 0) != -1 ) {
            numNodes++;
            for(int j = 0; j < NUM_NODES; j++) {
                if((amat(j, 0) >= 0) && (amat(j, i) == 1))
                numEdges ++;
                
            }
        }
    }
    numEdges /= 2;
    if(numEdges == (numNodes*(numNodes -1))/2) {
        return numNodes - 1;
    }
    
    for(int i = recDepth; i < 3; i++) {
        nodesChecked.push_back(nodesChecked[recDepth-1]);
    }
    int nodes[3];

    if(recDepth <= 3) {
        double val = data[nodesChecked[0]][nodesChecked[1]][nodesChecked[2]];
        if(val != -2) {
            numAvoided[numNodes]++;
            return val;
        }
        else {
            for(int i = 0; i < 3; i++) {
                nodes[i] = nodesChecked[i];
            }
            std::sort(nodes, nodes + 3);
        }
    }

    if(numNodes == 0) {
        retVal = -1;
    }
    
    else if (numNodes > 1) {
        if(numNodes == 3) {
            
        }
        for(int i = 0; i < NUM_NODES; i++) {
            if(amat(0, i) >= 0) {
                retVal += dimension(::unitSphere(amat, i), data, numAvoided);
            }
        }
    }
    
    if(numNodes > 1 ) {
        retVal = 1 + (retVal/numNodes);
    }
    if(recDepth <= 3) {
        do {
            
            data[nodes[0]][nodes[1]][nodes[2]] = retVal;
            
        } while(std::next_permutation(nodes, nodes+3));
        
    }
    
    
    return retVal;
    
    
    
}

void hGraph::oldCalcDimension() { //Calculates dimensionality recursively. See Knill, "On the dimensionanity and Euler Characteristic of Random Graphs"

    int kEdges = ((NUM_NODES)*(NUM_NODES - 1))/2;
    int kSum = 0;
    for(int i = 0; i < NUM_NODES; i++) {
        
        kSum += _degVector[i];
        
    }
    
    
    kSum /= 2;
    
    if(NUM_NODES == 0) {
        _dimension =  -1;
    }
    else if(kSum == kEdges) {
        
        //_dimension = NUM_NODES - 1;
        std::cout << NUM_NODES - 1 << std::endl;
        return;
    }
    else if (_numThreads == 1) {
        double dimen = 0.0;
        for(int i = 0; i < NUM_NODES; i++ ) {
            dimen += unitSphere(i).oldDimension(0, unitSphere(i).getSize(), false);
        }
        dimen = dimen/(static_cast<double>(NUM_NODES));
        _dimension = (dimen + 1);
        
    }
    
    else {
        int num;
        int extra;
        if(NUM_NODES < _numThreads) {
            num = 1;
            _numThreads = NUM_NODES;
            extra = 0;
        }
        else {
            num = NUM_NODES/_numThreads;
            extra = NUM_NODES%_numThreads; //extra nodes that don't divide evenly into number of threads.
            
        }
        int distributed = 0;
        std::vector<std::future<double>> futures;
        for(int i = 0; i < _numThreads -1; i++) {
            if(extra != 0) {
                distributed++;
                extra--;//If there are extra nodes, they are distributed to the various threads as they are called until there are no extra ones left.
                futures.push_back(std::async(std::launch::async, &hGraph::oldDimension, this, num*i + distributed-1, num*(i+1) + distributed, true) );
            }
            else {
                futures.push_back(std::async(std::launch::async, &hGraph::oldDimension, this, num*i + distributed, num*(i+1) + distributed, true) );

            }
        }
        double sum = oldDimension(num*(_numThreads-1) + distributed, NUM_NODES, true); //The main thread does some work as well.
        for(int i = 0; i <_numThreads -1; i++) {
            sum += futures[i].get();
        }
        double dimen = sum/(static_cast<double>(NUM_NODES));
        //_dimension = dimen + 1;
        std::cout << dimen + 1 << std::endl;
    }
    
    
    
}

double hGraph::oldDimension(int a, int b, bool multi) { //This function is NOT called by the user. It is used by the public function.
    int kEdges = ((NUM_NODES)*(NUM_NODES - 1))/2;   //Allows the function to be multithreaded with an arbitrary number of threads.
    int kSum = 0;
    for(int i = 0; i < NUM_NODES; i++) {
        
        kSum += _degVector[i];
        
    }
    kSum /= 2;
    
    if(NUM_NODES == 0) {
        return -1;
    }
    
    else if(kSum == kEdges) {
        return NUM_NODES -1;
    }
    
    else {
        double dimen = 0.0;
        for(int i = a; i < b; i++ ) {
            hGraph unit = unitSphere(i);
            int size = unit.getSize();
            dimen += unitSphere(i).oldDimension(0, size, false);
        }
        
        if(!multi) {
            dimen = dimen/(static_cast<double>(NUM_NODES));
            return (dimen + 1);
        }
        else {
            return dimen;
        }
        
    }
}

std::vector<int> hGraph::fractionalDimen() {
    std::vector<int> values = {0, 1};
    if(NUM_NODES == 0) {
        values[0] = -1;
        return values;
    }
    else {
        for(int i = 0; i < NUM_NODES; i++) {
            std::vector<int> temp = unitSphere(i).fractionalDimen();
            values = fractionAdd(values, temp);
        }
        values[1] *= NUM_NODES;
        values = fractionAdd(values, {1,1});
        return values;
        
    }
    
    
}

void hGraph::calcSpectralDimen() {
    Eigen::Matrix<unsigned long long int, Eigen::Dynamic, Eigen::Dynamic> walkMatrix;
    walkMatrix.resize(NUM_NODES, NUM_NODES);
    walkMatrix = _adjMatrix.cast <unsigned long long int> ();
    std::vector<double> temp;
    for(int i = 0; i < 500; i++) {
        double sum = 0;
        walkMatrix *= _adjMatrix.cast <unsigned long long int> ();
        for(int j = 0; j <  NUM_NODES; j++) {
            double partial = walkMatrix(j,j);
            double num = 0;
            for(int k = 0; k < NUM_NODES; k++) {
                num += walkMatrix(j,k);
            }
            partial /= num;
            sum += partial;
        }
        sum /= NUM_NODES;
        if(sum <= 0) {
            break;
        }
        temp.push_back(log(sum));
    }
    for(int i = 0; i < temp.size(); i++) {
        //double temp1 = temp[i] - temp[i-1];
        //double temp2 = log(i+2) - log(i+1);
        _spectralDimen.push_back(-2.0*(temp[i]/log(i+2)));
        
    }
    
}

void hGraph::calculateHausDimen() {
    MatrixXi temp;
    temp.resize(NUM_NODES, NUM_NODES);
    temp = getShortestPaths();
    int i = 2;
    while(true) {
        int numFound = 0;
        
        for (int j = 0; j < NUM_NODES; j++) {
            for (int k = 0; k < NUM_NODES; k++) {
                if(temp(j,k) == i) {
                    numFound += 1;
                }
            }
        }
        if(numFound == 0 ) {
            break;
        }
        _hausdorffDimen.push_back(log(numFound)/log(i));
        i++;
    }
    
    

    
    
}


std::vector<double> hGraph::autoGroupSize() {
    int m, n;
    graph g[MAXN*MAXM];
    int lab[MAXN], ptn[MAXN], orbits[MAXN];
    static DEFAULTOPTIONS_GRAPH(options);
    statsblk stats;
    n = NUM_NODES;
    m = SETWORDSNEEDED(n);
    EMPTYGRAPH(g, m, n);
    for(int i = 0; i < n; i++) {
        for(int j = i+1; j < n; j++) {
            if(isConnected(i,j)) {
                
                ADDONEEDGE(g, i, j, m);
            }
        }
    }
    densenauty(g, lab, ptn, orbits, &options, &stats, m, n, NULL);
    std::vector<double> temp;
    temp.push_back(stats.grpsize1);
    temp.push_back(stats.grpsize2);
    return temp;
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
    /*std::cout << "Number of cliques of each size" << std::endl;
    for(int i = 0; i < NUM_NODES; i++) {
        std::cout << i+1 << ": " << _numCliques[i] << std::endl;
        
    }*/
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

void hGraph::calculateEccen() {
    Eigen::Matrix<long unsigned int, Eigen::Dynamic, Eigen::Dynamic> walkMatrix;
    walkMatrix.resize(NUM_NODES, NUM_NODES);
    walkMatrix = _adjMatrix.cast<long unsigned int>();
    int length = 1;
    int numFound = 0;
    while(true) {
        for(int i = 0; i < NUM_NODES; i++) {
            
            if(_eccentricities[i] != 0) {
                continue;
            }
            
            for(int j = 0; j < NUM_NODES; j++ ) {
                if (walkMatrix(i,j) == 0 && (i != j)) {
                    break;
                }
                else if(j == (NUM_NODES-1)) {
                    _eccentricities[i] = length;
                    numFound++;
                }
            }
        }
        if((numFound == NUM_NODES) || (length == NUM_NODES)) {
            break;
        }
        walkMatrix *= _adjMatrix.cast<long unsigned int>();
        length++;
    }
    
}


/*int hGraph::pathL(int length, int a, int b) { //determines the number of paths of given length that connect vertices a and b. See arXiv:0801.0861 [hep-th], page 4
    int sum = 0;
    int prod = 0;
    for(int i = 0; i < NUM_NODES; i++) {
        prod = _adjMatrix(a,i) * _adjMatrix(i,b);
        sum += prod;
        
    }
    
    
    
    
}*/


//---------------------------END CALCULATIONS---------------------------//


//---------------------------I/O FUNCTIONS---------------------------//

void hGraph::toFile(std::ofstream &fs) const { //outputs hGraphs to a file in a CSV format. See documentation for information on formatting and usage.
    for(int i = 0; i < NUM_NODES; i++) {
        for(int j = 0; j < NUM_NODES; j++) {
            int temp = _adjMatrix(i,j);
            fs << temp << ",";
            
        }
        fs << "\n";
    }
    fs << "\n";
    for(int i = 0; i < NUM_NODES; i++) {
        fs << _degVector[i] << ",";
    }
    
    fs << "\n";
    
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
    _dimension = 0;
    _numCliques = std::vector<int> (size, 0);
    cliquesFound = false;
    
    
    for(int i = 0; i < NUM_NODES; i++) {
        for(int j = 0; j < size; j++) {
            _degVector[i] += _adjMatrix(i, j);
        }
    }
    _eccentricities.clear();
    for(int i = 0; i < size; i++) {
        _eccentricities.push_back(0);
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

////------------------------BEGIN INTPAIR CLASS IMPLEMENTATION------------------/

intPair::intPair(int pointA, int pointB) {
    a = pointA;
    b = pointB;
}

intPair::intPair() {
    a = -1;
    b = -1;
    
}

void intPair::print() {
    std::cout << "(" << a << ", " << b << ")";
    
}

std::vector<int> intPair::getPair() {
    std::vector<int> temp;
    temp.push_back(a);
    temp.push_back(b);
    return temp;
    
}

bool intPair::equals(intPair other) {
    std::vector<int> otherVector = other.getPair();
    if((otherVector[0] == a) && (otherVector[1] == b)) {
        return true;
    }
    else if((otherVector[0] == b) && (otherVector[1] == a)) {
        return true;
    }
    
    else
        return false;
    
    
}

/////-------------------------OTHER RESOURCE FUNCTIONS-------------------------/////

//---------------------------GRAPH GENERATORS---------------------------//

hGraph compGraph(int size) { //generates a complete graph of given size.
    
    MatrixXi adjMatrix = MatrixXi::Zero(size, size);
    for(int i = 0; i < size; i++) {
        for(int j = 0; j < size; j++) {
            if(j != i) {
                adjMatrix(i, j) = 1;
            }
        }
    }
    
    hGraph graph(size, adjMatrix);
    return graph;
    
    
    
    
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
    std::uniform_int_distribution<int> dist(max/4, (3*max)/4); //The graph will have at least 1/4th of the maximum number of possible edges
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

hGraph zeroGraph(int size) {
    MatrixXi adjMatrix = MatrixXi::Zero(size, size);
    hGraph graph(size, adjMatrix);
    return graph;

    
}

//---------------------------END GRAPH GENERATORS---------------------------//


//---------------------------I/O UTILITY FUNCTIONS---------------------------//



void readGraphFile(hGraph *** graphs, int &num) { //reads graphs from a CSV, returns a pointer to an array of hGraphs and sets variable num to the number of graphs read from the file.
                                   //See documentation for how to output graphs to CSV.
    int size;
    char cSize;
    std::string filename;
    std::ifstream input;
    std::cin.clear();
    while(true) {
        std::cout << "Input graph data filename: "; //gets input file name
        std::getline(std::cin, filename);
        input.open(filename);
        if(input.good()) {                  //makes sure file exists
            std::cout << "File found." << std::endl;
            break;
        }
        else {
            std::cout << "Invalid file. Check that the filename is spelled correctly and that it is in the main folder." << std::endl;
        }
    }
    
    
    std::string line;
    getline(input, line);
    size = stoi(line);
    int numRead = 0;

    
    int data[size*size];        //creates array that will temporarily hold adjacency matrix data;
    
    
    int numLines = size + 3;
    int count = 0;
    while(true) {
        getline(input,line); //gets the next line of the CSV file
        
        if(input.eof()) {  //checks if end of CSV file has been reached;
            break;
        }
        
        count++;
        if(count == numLines) {
            numRead++;
            count = 0;
        }
    }
    input.clear();
    input.seekg(0, std::ios::beg);
    getline(input, line);
    *graphs = new hGraph * [numRead];
    
    std::cout << "Reading in graph data..." << std::endl;
    for(int i = 0; i < numRead; i++) {
        MatrixXi adjMatrix = MatrixXi::Zero(size, size);
        for(int j = 0; j < numLines - 3; j ++) {
            getline(input,line); //gets the next line of the CSV file
            for(int k = 0; k < 2*size; k+= 2) {
                
                adjMatrix(j, k/2 ) = (line[k] - '0');
                
            }

            
        }
        
        (*graphs)[i] = new hGraph(size, adjMatrix);
        getline(input, line);
        getline(input, line);
        getline(input, line);
                
    }
    
    
    input.close();
    num = numRead;
    std::cout << "Number of graphs read: " << numRead << std::endl;
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

MatrixXi unitSphere(MatrixXi matrix, int node) {

    for (int i = 0; i < matrix.rows(); i++) {
        
        if(matrix(i, node) == 0) {      //Basic algorithm. Simply removes the row and column for a node from the adajecny matrix
            matrix(0, i) = -1;
            matrix(i, 0) = -1;
        }
        
    }

    matrix(node, 0) = -2;
    matrix(0, node) = -2;
    return matrix;
    
}

//---------------------------END MATRIX MUTATION FUNCTIONS---------------------------//

//---------------------------OTHER CALCULATION FUNCTIONS---------------------------//

bool isIsomorphic(hGraph graph1, hGraph graph2) {
    
    
    graph g1[MAXN*MAXM];
    graph g1b[MAXN*MAXM];
    graph g2[MAXN*MAXM];
    graph g2b[MAXN*MAXM];
    int lab[MAXN],ptn[MAXN],orbits[MAXN];
    static DEFAULTOPTIONS_GRAPH(options);
    statsblk stats;
    int n,m;
    options.writeautoms = FALSE;
    options.getcanon = TRUE;
    if(graph1.getSize() != graph2.getSize()) {
        return false;
    }

    n = graph1.getSize();
    m = SETWORDSNEEDED(n);
    EMPTYGRAPH(g1,m,n);
    EMPTYGRAPH(g1b,m,n);

    for (int i = 0; i < n; i++) {
        for(int j = i+1; j < n; j++) {
            if(graph1.isConnected(i,j)) {
                ADDONEEDGE(g1,i,j,m);
            }
        }
    }
    densenauty(g1,lab,ptn,orbits,&options,&stats,m,n,g1b);
    
    EMPTYGRAPH(g2,m,n);
    EMPTYGRAPH(g2b,m,n);
    for (int i = 0; i < n; i++) {
        for(int j = i+1; j < n; j++) {
            if(graph2.isConnected(i,j)) {
                ADDONEEDGE(g2,i,j,m);
            }
        }
    }
    densenauty(g2,lab,ptn,orbits,&options,&stats,m,n,g2b);
    
    int k;
    for(k = 0; k < m*(size_t)n; k++) {
        if(g1b[k] != g2b[k]) {
            break;
        }
    }
    
    if(k == m*(size_t)n) {
        return true;
    }
    else {
        return false;
    }
}

std::vector<int> fractionAdd(std::vector<int> fractA, std::vector<int> fractB) {
    
    int LCD = fractA[1]*fractB[1]/gcd(fractA[1], fractB[1]);
    int NUM = fractA[0]*(LCD/fractA[1]) + fractB[0]*(LCD/fractB[1]);
    std::vector<int> result = {NUM, LCD};
    return result;
    
}

void fractReduce(std::vector<int> &frac) {
    if(frac.size() != 2) {
        std::cout << "Error: vector passed is not a fraction " << std::endl;
        exit(1);
    }
    int GCD = gcd(abs(frac[0]), frac[1]);
    frac[0] /= GCD;
    frac[1] /= GCD;
    
}


int gcd(int numA, int numB) {
    int num1;
    int num2;
    if(numA > numB) {
        num1 = numA;
        num2 = numB;
    }
    else {
        num1 = numB;
        num2 = numA;
    }
    while((num1 % num2) != 0) {
        int temp = num1 % num2;
        num1 = num2;
        num2 = temp;
        
    }

    return num2;
}



//---------------------------END OTHER CALCULATION FUNCTIONS---------------------------//

/////-------------------------END OTHER RESOURCE FUNCTIONS-------------------------/////


