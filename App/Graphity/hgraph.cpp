//
//  hGraph.cpp
//
//
//  Created by Patrick on 10/14/17.
//
//

#include <algorithm>
#include "hGraph.h"
//#include "absHamiltonian.h"
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

hGraph::hGraph(int size):NUM_NODES(size) { //Basic constructor. Creates an graph of given size with no connections. Most often used of the constructors
 _adjMatrix = MatrixXi::Zero(size, size);
 _degVector = std::vector<int>(size , 0);
 _numCliques = std::vector <int> (size, 0);
 _eccentricities = std::vector<int> (size, 0); //Initializes this vector with zeroes


}

hGraph::hGraph() {  //default constructor. Creates placeholder. Should not be used often.

 _adjMatrix = MatrixXi::Zero(1, 1);
 _degVector = std::vector<int> {0};
 _hamiltonian = 0.0;

}

hGraph::hGraph(int size, MatrixXi adjMatrix):NUM_NODES(size) { //Takes a size and a previously defined adjacency matrix and creates an hGraph object.
 _adjMatrix = MatrixXi::Zero(size, size); //Matrix has to be created and assigned in seperate steps
 _adjMatrix = adjMatrix;
 _degVector = std::vector<int> (size, 0);
 _hamiltonian = 0.0;
 _dimension = 0; //this is not initialized by default due to the high computing cost.
 _numCliques = std::vector<int> (size, 0);
 _eccentricities = std::vector<int> (size, 0);
 for(int i = 0; i < NUM_NODES; i++) {
     for(int j = 0; j < NUM_NODES; j++) {
         _degVector[i] += _adjMatrix(i, j);
     }
 }



}

//---------------------------END CONSTRUCTORS---------------------------//


//---------------------------DESTRUCTOR---------------------------//


hGraph::~hGraph() { //Since the hGraph class does not include any pointers, there is no work for this destructor to do.

}

//---------------------------END DESTRUCTOR---------------------------//

//---------------------------GETTERS---------------------------//

double hGraph::getDimension() { //simple getter function. See below for the dimensionality algorithm.
 if(!dimensionFound) { //Checks if the dimenison has already been calculated
     calcDimension();
     dimensionFound = true;
 }
 return _dimension;

}

std::vector<double> hGraph::getSpectralDimen() { //simple getter function. See below for the spectral dimensionality algorithm
 if (_spectralDimen.size() == 0) {
     calcSpectralDimen();
 }
 return _spectralDimen;

}

double hGraph::getHam() const {  //simple getter function
 return _hamiltonian;
}

int hGraph::getEulerChar() { //simple getter function. See below for the euler characteristic agorithm.
 calcEulerChar();
 return _eulerChar;
}

int hGraph:: getSize() const{ //simple getter function
 return NUM_NODES;
}

std::vector<int> hGraph::numCliques() { //simply outputs the number of cliques (NOT maximal cliques) of any given size. The cliques are stored in an vector attribute called "
 countCliques(); //this function checks if the cliques have already been found, and simply terminates if so.
 return _numCliques;

}

int hGraph::getDegree(int node) const { //simply gets the degree of a given node.
 if(node < 0 || node >= NUM_NODES) { //Recall that C++ is zero indexed, the nodes of a graph of size 10 will be labled 0...9
     std::cerr << "Critical error, only acceptable node values are between 0 and " << NUM_NODES - 1 << std::endl;
     exit(4);
 }

 return _degVector[node];

}

bool hGraph::isConnected(int row, int column) const { //checks if two nodes are connected.
 if(row < 0 || row >= NUM_NODES || column < 0 || column > NUM_NODES) { //Recall that C++ is zero indexed, the nodes of a graph of size 10 will be labled 0...9
     std::cerr << "Critical error, only acceptable node values are between 0 and " << NUM_NODES - 1 << std::endl;
     exit(4);
 }


 if(_adjMatrix(row, column) == 1) {
     return true;
 }

 return false;

}

std::vector<int> hGraph::getEccentricity() { //eccentricity of a node is how minimum number of walks that can get from one node to any other node
 if(_eccentricities[0] == 0 && _eccentricities[NUM_NODES/2] == 0) {
     calculateEccen();
 }
 return _eccentricities;

}

int hGraph::getDiameter() {
 std::vector<int> temp = getEccentricity(); //The dimater of the graph is the max ecentricity of one of its nodes.
 int val = 0;
 for(int i = 0; i < temp.size(); i++ ) {
     if(temp[i] > val) {
         val = temp[i];
     }
 }

 return val;

}

std::vector<double> hGraph::getHausdorffDimen() { //getter function. See below for calculation of the hausdorff dimension

 if(_hausdorffDimen.size() == 0) {
     calculateHausDimen();

 }

 return _hausdorffDimen;

}



//---------------------------END GETTERS---------------------------//

//---------------------------SETTERS---------------------------//

void hGraph::setHamiltonian(double val) { //sets the value of the hamiltonian.
 _hamiltonian = val;
}

void hGraph::setThreads(int threads) { //Note, this does NOT check to ensure that the given number of threads can actually be used. That is done by multithreaded algorithms.
 _numThreads = threads;
}

void hGraph::flipEdge(int nodeA, int nodeB) {
 if(isConnected(nodeA, nodeB)) { //Turns on an edge if it doesn't exist, or turns it off if it does. Updates node degress accordingly.
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

 cliquesFound = false; //In general, if an edge if flipped this information is lost.
 dimensionFound = false; //There is no closed-form algorithm (that I am aware) that determines the new dimensionality if a single edge is removed.
 curvatureAtNode = std::vector<double>(NUM_NODES, 0);


}

void hGraph::flipEdge(std::vector<int> nodeA, std::vector<int> nodeB) { //Multiple edges can be flipped at once. The two vectors correspond to pairs of nodes.
 for(int i = 0; i < nodeA.size(); i++) {                             //This is used primarily in the monte-carlo algorithm.
     flipEdge(nodeA[i], nodeB[i]);
 }
}

void hGraph::acceptPartial(double partial) { //Adds the value passed to the function to the value of the hamiltonian. Used primarily in the Monte Carlo simulation
 _hamiltonian += partial;

}

//---------------------------END SETTERS---------------------------//


//---------------------------CALCULATIONS---------------------------//

void hGraph::calcEulerChar() { //Knill presents two equivalent definitions of the Euler Characteristic of a graph.
 double sum = 0;            //Definition used here is the sum of the curvatures at every node, because it is easily parallelizable
 int numFound = 0;
 for(int i = 0; i < NUM_NODES; i++) {
     if(curvatureAtNode.size() == i) {
         curvatureAtNode.push_back(0);
     }
     else if(curvatureAtNode[i] != 0) {
         numFound++; //Counts how many curvatures have already been found.
     }


 }
 if(_numThreads == 1) {
     for(int i = 0; i < NUM_NODES; i++) {
         if(curvatureAtNode[i] == 0) { //Checks if the curvature at that node has already been calculated.
             curvatureAt(i); //While curvatureAt does return the value, it also stores it in the vector curvatureAtNode,
         }                   //which is an attribute of the hGraph class. The sum is calculated at the bottom of the function
     }
 }


 else if(_numThreads >= NUM_NODES) { //If there are as many or more threads than their are nodes.

     std::vector<std::future<double>> futures;
     bool mainThread = (numFound == 0 ? true : false); //If some of the curvatures have already been found, the main thread will not need to do a calculation.
     int i = 0;                                        //If you are unfamiliar with the syntax in the above line, look up "c++ ternary operator"
     while(true) {
         if(i == (NUM_NODES - 1)) {
             if(mainThread) {
                 break;
             }
             else {
                 futures.push_back(std::async(std::launch::async, &hGraph::curvatureAt, this, i));
                 break;
             }
         }
         if(curvatureAtNode[i] == 0) {
             futures.push_back(std::async(std::launch::async, &hGraph::curvatureAt, this, i));

         }
         i++;
     }

     double tempA = 0;
     double count = 0;

     if(mainThread) {
         curvatureAt(NUM_NODES - 1);
     }
     for(int i = 0; i < futures.size(); i++) {
             futures[i].get(); //Again. While this will return a value, it is already stored elsewhere so it is unnecessary keep it here.
     }
 }

 else { //This is what runs if there are more nodes than their are threads.
     int numLaunches = (NUM_NODES-numFound)/_numThreads;//Number of times threads will have to be launched (i.e. for 30 nodes and 10 threads each thread will process 3 nodes)
     int extra = (NUM_NODES-numFound) % _numThreads; //Nodes not covered by the aboce (i.e. for 32 nodes and 10 threads, there will be two leftover nodes).
     numLaunches += (extra == 0 ? 0 : 1); //There has to be one more launch if there are extra nodes.
     int found = 0; //Nodes are processed in numerical order. This variable tracks how many nodes have been processed that had ALREADY been calculated.
     for(int i = 0; i < numLaunches; i++) {
         std::vector<std::future<double>> futures; //Holds the thread objects

         if(i != (numLaunches -1) || extra == 0 ) { //If this particulary "launch" will use all available threads
             int j = 0;
             while(true) {
                 int index = i*_numThreads + j + found;
                 if(curvatureAtNode[index] == 0) { //if the value has not yet been found
                     futures.push_back(std::async(std::launch::async, &hGraph::curvatureAt, this, index));
                     j++;
                 }
                 else {
                     found++;
                 }

                 if(j == _numThreads - 1) { //Allowx the main thread to do some of the work.
                     break;
                 }

             }
             curvatureAt((i+1)*_numThreads + found -1); //main thread does work.

         }
         else { //Calculates for the "extra" nodes.
             int j = 0;
             while(true) {
                 int index = i*_numThreads + j + found;
                 if(curvatureAtNode[index] == 0) {
                     futures.push_back(std::async(std::launch::async, &hGraph::curvatureAt, this, index));
                     j++;
                 }
                 else {
                     found++;
                 }

                 if(j == extra) {
                     break;
                 }

             }


         }
         for(int j = 0; j < futures.size(); j++) {
             futures[j].get();

         }

     }


 }
 for(int i = 0; i < NUM_NODES; i++) {
     sum += curvatureAtNode[i]; //Calculates the sum from the vector that was populated previously.

 }

 _eulerChar = round(sum);
 //The Euler Characteristic is always an integer, but the individual curvatures may not be. As such, it is possible that
 //a rounding error will cause the sum to be very very close to, but not quite equal to, the correct value.
 //Assigning a double to an int truncates the decimal, so rounding the value ensures it is correct.

}

double hGraph::curvatureAt(int node) {

 if(curvatureAtNode.size() != 0 && curvatureAtNode[node] != 0) {
     return curvatureAtNode[node];
 }

 else if(curvatureAtNode.size() == 0) {

     curvatureAtNode = std::vector<double>(NUM_NODES, 0);

 }

 hGraph temp = unitSphere(node);
 if(temp.getSize() == 0) {
     curvatureAtNode[node] = 1;
     return 1;
 }
 std::vector <int> tempVector = temp.numCliques();
 double sum = 0;
 for(int i = 0; i < temp.getSize(); i++) {
     double val = double(tempVector[i])/(i+2);
     if(i%2 == 0) {
         sum -= val;
     }
     else {
         sum += val;

     }
 }

 curvatureAtNode[node] = sum+1;
 return sum + 1;


}


MatrixXi hGraph::getShortestPaths() {
 MatrixXi pathMatrix(NUM_NODES, NUM_NODES); //uses the Floyd-Warshall algorithm to determine the shortest path between all nodes in a graph.
 pathMatrix = _adjMatrix;                   //The entry i, j  of the resultant matrix will be the minimum number of walks from i to j
 for(int i = 0; i < NUM_NODES; i++) {
     for (int j = 0; j < NUM_NODES; j++ ) {

         if((i != j) && (pathMatrix(i,j) == 0)) {
             pathMatrix(i,j) = 10000000; //A stand in for infinity. Two nodes that CANNOT be walked between will have this value (unless the nodes are the same node)
         }                               //Both Eigen and C++ have infinity as a value, which was causing issues, so I just used a big number.

     }
 }


 for(int k = 0; k < NUM_NODES; k++) {    //Calculates the actual values
     for(int i = 0; i < NUM_NODES; i++) {
         for(int j = 0; j < NUM_NODES; j++) {

             if(pathMatrix(i,j) > (pathMatrix(i,k) + pathMatrix(j, k))) {
                 pathMatrix(i,j) = (pathMatrix(i,k) + pathMatrix(j,k));
             }


         }
     }
 }


 for(int i = 0; i < NUM_NODES; i++) {
     for(int j = 0; j < NUM_NODES; j++) {
         if(pathMatrix(i, j) == 10000000) {
             pathMatrix(i, j) = -1;
         }
     }
 }

 return pathMatrix;


}

hGraph hGraph::compliment() const { //Returns an hGraph object where the adjacency matrix is the compliment of the original hGraph's adjacency matrix.
 hGraph temp(NUM_NODES, _adjMatrix);
 for(int i = 0; i < NUM_NODES; i++) {
     for(int j = i+1; j < NUM_NODES; j++) { //Using nested loops where the inner loop
         temp.flipEdge(i, j);
     }
 }

 return temp;


}

double hGraph::getAvgDegree() const { //Returns the average node degree
 double avg = 0;
 for(int i = 0; i < NUM_NODES; i++) {
     avg += _degVector[i];
 }
 return avg/NUM_NODES;
}


void hGraph::calcDimension() { //The base algorithm can be found in "On the Dimensionality and Euler Characteristic of Random Graphs" by O.Knill

 int num = 0;               //A more detailed explanation of this particular implementation can be found in the users guide
 int maxEdges = NUM_NODES*(NUM_NODES-1) / 2;
 int edges = 0;
 for(int i = 0; i < NUM_NODES; i++) {
     edges+= getDegree(i);
 }
 edges /= 2;
 if(edges == maxEdges) {
     _dimension = NUM_NODES -1;
     return;
 }


 std::vector<int> * trip;
 trip = new std::vector<int>[3];
 double *** values;
 values = new double ** [NUM_NODES];
 int triples = 0;


 for(int i = 0; i < NUM_NODES; i++) {

     values[i] = new double * [NUM_NODES-i];

     for(int j = 0; j < NUM_NODES - i; j++) {

         values[i][j] = new double[NUM_NODES-j];
         for(int k = 0; k < NUM_NODES - j; k++) {

             values[i][j][k] = -2.0;
             if(i < i + j && i + j < j + k) {
                 if(isConnected(i, i+j) && isConnected(i + j, j + k)) { //ensures the depth-3 sphere in question actually exists
                     trip[0].push_back(i);
                     trip[1].push_back(i + j);
                     trip[2].push_back(j + k);
                     triples++;
                 }

             }
         }
     }
 }
 if(_numThreads == 1) {
     dimension(_adjMatrix, values, trip, 0, trip[0].size(), true); //function calculates the values of the depth-3 sphere
 }

 else {


     int tripSize = trip[0].size();
     if(_numThreads > tripSize) { //Ensures that threads are not created if there's no work for them to do.
         _numThreads = tripSize;
     }

     int extra = tripSize%_numThreads;
     int num = tripSize/_numThreads;
     int distributed = 0;
     std::vector<std::future<double>> futures;

     for(int i = 0; i < _numThreads -1; i++) { //Distributes the depth-3 spheres calculations across all of the available threads.
                                               //If the number of depth-3 spheres does not divide evenly into the number of threads, the extras are distributed
         if(extra > 0) {
             distributed++;
             extra--;
             int low = i*num + distributed - 1;
             int high = (i+1)*num + distributed;
             futures.push_back(std::async(std::launch::async, &hGraph::dimension, this, _adjMatrix, values, trip, low, high, true));
         }

         else {
             int low = num*i + distributed;
             int high = num*(i+1) + distributed;

             futures.push_back(std::async(std::launch::async, &hGraph::dimension, this, _adjMatrix, values, trip, low, high, true));

         }


     }
     dimension(_adjMatrix, values, trip, num*(_numThreads-1) + distributed, tripSize, true); //the main thread also does osme work
     for(int i = 0; i < _numThreads -1; i++) {
         futures[i].get();
     }

 }
 delete [] trip; //The triples are no longer actually needed, so the pointer is deleted

 //This value will hold the sum of all the spheres of depth 1, and will be updated as they are calculated.
 for(int i = 0; i < NUM_NODES; i++) { //This loop calculates all the depth-2 spheres that will go into the individual depth-1 spheres

     for(int j = i+1; j < NUM_NODES; j++) {


         double sph2dimen = 0;
         if(i == j || !isConnected(i, j) ) {
             continue;
         }
         double numNodes = 0;


         for(int k = 0; k < NUM_NODES; k++) {
             if(isConnected(i,k) && isConnected(j,k)) {
                 numNodes++;
             }
         }

         if(numNodes == 0) {
             sph2dimen = -1;
         }
         else if(numNodes == 1) {

         }

         else {
             double sumA = 0;
             for(int k = 0; k < NUM_NODES; k++) {
                 if(i == k || j == k ||  !isConnected(i, k) || !isConnected(j,k)) {
                     continue;
                 }
                 int spheres[3] = {i, j, k};
                 std::sort(spheres, spheres + 3);
                 double val = values[spheres[0]][spheres[1] - spheres[0]][spheres[2] - spheres[1]];
                 sumA += values[spheres[0]][spheres[1] - spheres[0]][spheres[2] - spheres[1]];
             }
             sph2dimen = 1 + (sumA/numNodes);

         }
         values[i][j-i][0] = sph2dimen;

     }

 }

 double sphSum1 = 0;
 for(int i = 0; i < NUM_NODES; i++) { //calculates the 1-spheres that will go into the finan dimensionality calculation
     double sphSum2 = 0;
     int numNodes = 0;

     for(int j = 0; j < NUM_NODES; j++) {
         if(!isConnected(i, j)) {
             continue;
         }
         else {
             numNodes++;
             int tempPair[2] = {i , j};
             std::sort(tempPair, tempPair + 2);
             double val = values[tempPair[0]][tempPair[1] - tempPair[0]][0];
             sphSum2 += values[tempPair[0]][tempPair[1] - tempPair[0]][0];
         }
     }
     if(numNodes == 0) {
         sphSum1 += -1;
     }
     else if (numNodes == 1) {

     }
     else {
         double dimen1 = 1 + sphSum2/numNodes;
         sphSum1 += 1+ (sphSum2/numNodes);
     }

 }
 _dimension = 1 + sphSum1/NUM_NODES; //final result


 for(int i = 0; i < NUM_NODES; i++) { //deletes pointers


     for(int j = 0; j < NUM_NODES - i; j++) {

         delete [] values[i][j];
     }

     delete [] values[i];

 }
 delete [] values;

}

double hGraph::dimension(MatrixXi amat, double *** data, std::vector<int> * trp, int lowerBound, int upperBound, bool init)  { //utility function used by main calcDimension function
 if(init) { //calculates all 3-sphres needed
     MatrixXi tempAmat(NUM_NODES, NUM_NODES);

     for(int i = lowerBound; i < upperBound; i++) {
         tempAmat = ::unitSphere(amat, trp[0][i]);
         tempAmat = ::unitSphere(tempAmat, trp[1][i]);
         tempAmat = ::unitSphere(tempAmat, trp[2][i]);
         double val = dimension(tempAmat, data, trp, 0, 0, false);
         data[trp[0][i]][trp[1][i] - trp[0][i]][trp[2][i] - trp[1][i]] = val;


     }

     return -2;

 }


 else { //This is uesd for recursive calls below the 3-sphree level
     int numNodes = 0;
     int numEdges = 0;
     for(int i = 0; i < NUM_NODES; i++) {
         if(amat(0,i) >= 0) {
             numNodes++;
             for(int j = 0; j < NUM_NODES; j++) {
                 if((amat(0,j) >= 0) && (amat(i,j) == 1)) {
                     numEdges++;
                 }
             }
         }
     }
     numEdges /= 2;
     if(numNodes == 0) {
         return - 1;
     }

     else if(numNodes == 1 || numEdges == 0) {
         return 0;
     }
     else if(numEdges == (numNodes*(numNodes -1))/2) {

         return numNodes -1;


     }

     else {
         double sum = 0;
         for (int i = 0; i < NUM_NODES; i++) {
             if(amat(0, i) >= 0) {
                 sum += dimension(::unitSphere(amat,i), data, trp, 0, 0, false);
             }
         }
         return 1 + (sum/numNodes);
     }
 }
}

void hGraph::oldCalcDimension() { //Calculates dimensionality recursively. See Knill, "On the dimensionanity and Euler Characteristic of Random Graphs"
                               //This algorithm is no longer in used, but is left here in case it comes up

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

     _dimension = NUM_NODES - 1;
     return;
 }
 else if (_numThreads == 1) {
     double dimen = 0.0;
     for(int i = 0; i < NUM_NODES; i++ ) {
         dimen += unitSphere(i).oldDimension(0, unitSphere(i).getSize(), false);
     }
     dimen = dimen/(static_cast<double>(NUM_NODES));
     _dimension = 1 + dimen;

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
     _dimension = dimen + 1;
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

std::vector<int> hGraph::getFractionalDimen() { //Single-threaded function that outputs result as a fraction. The first item in the vector is the numerator, second is denomenator
 std::vector<int> values = {0, 1};
 if(NUM_NODES == 0) {
     values[0] = -1;
     return values;
 }
 else {
     for(int i = 0; i < NUM_NODES; i++) {
         std::vector<int> temp = unitSphere(i).getFractionalDimen();
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
 walkMatrix = _adjMatrix.cast <unsigned long long int> (); //spectral dimension creates very large values in the matrix
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
     if(sum <= 0) { //Indicates that the numbers have become too big
         break;
     }
     temp.push_back(log(sum));
 }
 for(int i = 0; i < temp.size(); i++) {
     //double temp1 = temp[i] - temp[i-1];
     //double temp2 = log(i+2) - log(i+1);
     _spectralDimen.push_back(-2.0*(temp[i]/log(i+2))); //this is a vector for different walk lengths. In theory, the spectral dimension is the value when i -> infinity

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


std::vector<double> hGraph::autoGroupSize() const {
 //Calculates the size of the automorphism group using the Nauty/Traces library. See user guide.

 int m, n;
 graph g[MAXN*MAXM];
 int lab[MAXN], ptn[MAXN], orbits[MAXN];
 static DEFAULTOPTIONS_GRAPH(options);
 statsblk stats;
 n = NUM_NODES;
 m = SETWORDSNEEDED(n);
 EMPTYGRAPH(g, m, n); //Container for use in Nauty
 for(int i = 0; i < n; i++) {
     for(int j = i+1; j < n; j++) {
         if(isConnected(i,j)) { //copies the graph into Nauty's container.

             ADDONEEDGE(g, i, j, m);
         }
     }
 }
 densenauty(g, lab, ptn, orbits, &options, &stats, m, n, NULL);
 std::vector<double> temp;
 temp.push_back(stats.grpsize1); //Returns the value of the automorphism group size as a vector in scientific notation
 temp.push_back(stats.grpsize2); //First item is the base, second item is the exponent of the 10
 return temp;
}


hGraph hGraph::unitSphere(int node) const { //outputs the unit sphere around a given. Node the united sphere contains
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

void hGraph::countCliques() { //Declares necessary objects for initial call of cliqueSearch, then counts all cliques. Uses a moified version of the Bron-Kerbosch algorithm
 if(cliquesFound) {
     return;
 }
 std::vector<int> vectorR(0);
 std::vector<int> vectorP(NUM_NODES);
 for(int i = 0; i < NUM_NODES; i++) {
     vectorP[i] = i;
 }

 cliqueSearch(vectorR, vectorP);

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

void hGraph::calculateEccen() { //Calculates the eccentricity (distances to farthest node) for every node in the graph.
 Eigen::Matrix<long unsigned int, Eigen::Dynamic, Eigen::Dynamic> walkMatrix;
 walkMatrix.resize(NUM_NODES, NUM_NODES);
 walkMatrix = _adjMatrix.cast<long unsigned int>();
 int length = 1;
 int numFound = 0;
 while(true) { //Simply raises the matrix to a new power and then checks if each node has an entry for every other node
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


//---------------------------END CALCULATIONS---------------------------//


//---------------------------I/O FUNCTIONS---------------------------//

void hGraph::toFile(std::ofstream &fs) const { //outputs hGraphs to a file in a CSV format. See documentation for information on formatting and usage.
 for(int i = 0; i < NUM_NODES; i++) {
     for(int j = 0; j < NUM_NODES; j++) { //outputs adjacency matrix
         int temp = _adjMatrix(i,j);
         fs << temp << ",";

     }
     fs << "\n";
 }
 fs << "\n";
 for(int i = 0; i < NUM_NODES; i++) {
     fs << _degVector[i] << ","; //outputs node degrees
 }

 fs << "\n";

 fs << _eulerChar << "," << _hamiltonian << std::endl; //outputs hamiltonian
}

std::ofstream &operator << (std::ofstream &fs, const hGraph &rhs)  { //allows use of << operator to output hGraphs to files.
 rhs.toFile(fs);
 return fs;
}


void hGraph::toStream(std::ostream &os) const { //outputs hGraphs in a human-readable format for standard output. See documentation for information on formatting and usage.

 os << "Adjacency Matrix: " << std::endl << _adjMatrix << std::endl;
 os << "Node Degrees: " << std::endl;
 for(int i = 0; i < NUM_NODES; i++) {
     os << _degVector[i] << std::endl;
 }
 os << std::endl;
 os << "Dimensionality: " << _dimension << std::endl;
 os << "Euler Characteristic: " << _eulerChar << std::endl;
 os << "Value of most recently used hamiltonian: " << _hamiltonian << std::endl;

}

std::ostream &operator << (std::ostream &os, const hGraph &rhs)  { //Allows use of << operator to output hGraphs to standard output
 rhs.toStream(os);
 return os;
}



//---------------------------END I/O FUNCTIONS---------------------------//


//---------------------------HGRAPH RESOURCE FUNCTIONS---------------------------//



void hGraph::setMatrix(int size, MatrixXi data) { //resets the matrix based on new data. Should not be neccesary to use often.
 NUM_NODES = size;
 _adjMatrix.resize(size, size);
 _adjMatrix = data;
 _degVector = std::vector<int>(size, 0);
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

bool hGraph::operator != (hGraph & rhs) const {

 return !(*this == rhs);

}


bool hGraph::operator == (hGraph & rhs) const {

 if(this->getSize() != rhs.getSize()) {
     return false;
 }

 return isIsomorphic(*this, rhs);

}

//---------------------------END HGRAPH RESOURCE FUNCTIONS---------------------------//

/////-------------------------END HGRAPH CLASS IMPLEMENTATION-------------------------/////

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

hGraph randomGraph(int size, double fillFrac) { //generates a random graph of a given size. First argument gives number of nodes. Second gives filling fraction
                                                //If the second argument is zero, the filling fraction is also determined randomly.

    unsigned int max = size*(size-1)/2; //maximum number of edges

    std::random_device rd; //c++11 random number generator
    std::mt19937 gen(rd());



    int fill[max];
    for(int i = 0; i < max; i++ ) { //creates an array that can be mapped to an adjacency matrix and sets all items to zero.
        fill[i] = 0;
    }
    unsigned int fillNum;
    if(fillFrac == 0) { //Initialize filling fraction randomly if not provided by user.
        std::uniform_int_distribution<int> dist(max/4, (3*max)/4); //The graph will have at least 1/4th of the maximum number of possible edges, and at most 3/4ths
        fillNum = dist(gen);                                 //Ensures graphs are dense enough to be useful.

    }
    else {
        fillNum = ceil(fillFrac*max);
    }

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

hGraph zeroGraph(int size) { //generates a graph with no edges.
    MatrixXi adjMatrix = MatrixXi::Zero(size, size);
    hGraph graph(size, adjMatrix);
    return graph;


}

//---------------------------END GRAPH GENERATORS---------------------------//


//---------------------------I/O UTILITY FUNCTIONS---------------------------//

//the triple pointer is a C++ requirement regarding memory management in external functions

void readGraphFile(hGraph *** graphs, int &num) { //Takes in the address to a double (NOT TRIPLE) pointer to an hGraph, as well as an integer variable.
    int size;                                     //After function runs, graphs will be an array of pointers to graphs
    char cSize;                                   //variable num will contain the number of graphs read.
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

MatrixXi unitSphere(MatrixXi matrix, int node)  {

    //This function is used by the newer version of the dimensionality algorithm.
    //a -1 in the header indicates the node is not in the sphere, a -2 denotes the node the sphere was taken around.

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

bool isIsomorphic(hGraph graph1, hGraph graph2) { //Takes in two graphs and determines if they are isomorphic
                                                  //using the Nauty/Traces library
    graph g1[MAXN*MAXM];                          //Spcific logic is not important
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

//These last three functions were used for testing purposes. Does operations on fractions passed as vectors.

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



