   //
//   hGraph.h
//   MULTITHREADED VERSION
//
//  Created by Patrick Wells on 10/14/17.
//
//

#define TIMING

#ifdef TIMING
#define INIT_TIMER auto start = std::chrono::high_resolution_clock::now();
#define START_TIMER  start = std::chrono::high_resolution_clock::now();
#define STOP_TIMER(name)  std::cout << "RUNTIME of " << name << ": " << \
`).count() << " ms " << std::endl;
#else
#define INIT_TIMER
#define START_TIMER
#define STOP_TIMER(name)
#endif

#ifndef _hGraph_h
#define _hGraph_h

#include <iostream>
#include <fstream>
#include <vector>
#include <random>
#include <thread>
#include <future>
#include <limits>
#include <cmath>


#define MAXN 500
#include "nauty.h"


#include "../lib/eigen/Eigen/Dense"

using Eigen::MatrixXi;
class absHamiltonian;

class hList;
class hNode;

class hGraph {
    
private:
    int _numThreads = 1;
    int NUM_NODES;
    int _eulerChar = 0;
    double _dimension;
    std::vector<double> _spectralDimen;
    std::vector<double> _hausdorffDimen;
    double _hamiltonian;
    bool cliquesFound = false;
    
    std::vector<int> _eccentricities;


    MatrixXi _adjMatrix;
    Eigen::VectorXi _degVector;
    std::vector <int> _numCliques;
    
    void toStream(std::ostream &os) const;
    void toFile(std::ofstream &fs) const;
    double oldDimension(int a, int b, bool multi);     //Single threaded version for recursive calls
    double dimension(MatrixXi amat, double *** data, std::vector<int> * trp, int lowerBound, int upperBound, bool init);
    void calculateEccen();
    void calculateHausDimen();



    
public:
    hGraph(int size, MatrixXi adjMatrix);
    hGraph(int size);
    hGraph();
    ~hGraph();
    hGraph compliment();
    void print(); //depricated
    void setMatrix(int size, MatrixXi data);
    void setHamiltonian(double val);
    hGraph unitSphere(int node);
    void calcDimension();
    void oldCalcDimension(); //Multithreaded version
    std::vector<int> fractionalDimen();
    void calcSpectralDimen();
    void calcEulerChar();
    double avgDegree();
    void accept(absHamiltonian &ham);
    void flipEdge(int nodeA, int nodeB);
    void flipEdge(std::vector<int> nodeA, std::vector<int> nodeB);
    void acceptPartial(double partial);
    void setThreads(int threads);
    std::vector<int> getEccentricity();
    int getDiameter();
    std::vector<double> autoGroupSize();
    MatrixXi _paths;


    void numCliques();
    MatrixXi getShortestPaths();

    double getDimension();
    std::vector<double> getSpectralDimen();
    std::vector<double> getHausdorffDimen();
    
    int getDegree(int node);
    int getSize();
    int getEulerChar();
    bool isConnected(int row, int column);
    void countCliques();
    void cliqueSearch(std::vector<int> R, std::vector<int> P);
    double getHam();
    friend std::ostream &operator << (std::ostream &os, const hGraph &rhs);
    friend std::ofstream &operator <<(std::ofstream &fs, const hGraph &rhs);
};

//Forward declaration of hList class. Used to sort graphs based on their hamiltonian.

class hList {
private:
    hNode* _head;
    double _energy;
public:
    hList();
    ~hList();
    bool isEmpty();
    void clear();
    void prepend(hGraph * graph);
    void print();
    double getEnergy();
    void toFile(std::ofstream &fs) const;
    friend std::ostream &operator <<(std::ostream &os, const hList &rhs);
    friend std::ofstream &operator <<(std::ofstream &fs, const hList &rhs);

};

//forward declaration of hNode classe, used in hList

class hNode {
private:
    friend class hList;
    hNode * _next;
    hGraph * _data;
public:
    hNode(hGraph * graph);
    ~hNode();
    void clear();
    void print();
    void toFile(std::ofstream &fs) const;
};

class intPair {
private:
    int a;
    int b;
public:
    intPair(int pointA, int pointB);
    intPair();
    std::vector<int> getPair();
    
    void print();
    bool equals(intPair other);
};


//Forward declaration of needed resource functions

void readGraphFile(hGraph *** graphs, int &num);

MatrixXi unitSphere(MatrixXi matrix, int node);
void removeColumn(Eigen::MatrixXi& matrix, unsigned int colToRemove);
void removeRow(Eigen::MatrixXi& matrix, unsigned int rowToRemove);
bool isIsomorphic(hGraph graph1, hGraph graph2);
hGraph randomGraph(int size, double fillFrac);
hGraph compGraph(int size);
hGraph zeroGraph(int size);
std::vector<int> fractionAdd(std::vector<int> fractA, std::vector<int> fractB);
void fractReduce(std::vector<int> &frac);
int gcd(int numA, int numB);




#endif /* _hGraph_h */
