  //
//   hGraph.h
//   MULTITHREADED VERSION
//
//  Created by Patrick Wells on 10/14/17.
//
//

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
    double dimension(int a, int b, bool multi);     //Single threaded version for recursive calls
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
    void calcDimension(); //Multithreaded version
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


//Forward declaration of needed resource functions

void readGraphFile(hGraph *** graphs, int &num);

void removeColumn(Eigen::MatrixXi& matrix, unsigned int colToRemove);
void removeRow(Eigen::MatrixXi& matrix, unsigned int rowToRemove);
bool isIsomorphic(hGraph graph1, hGraph graph2);
hGraph randomGraph(int size);
hGraph compGraph(int size);
hGraph zeroGraph(int size);





#endif /* _hGraph_h */
