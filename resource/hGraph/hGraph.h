//
//   hGraph.h
//  
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

#include "../lib/eigen/Eigen/Dense"

using Eigen::MatrixXi;
class absHamiltonian;

class hList;
class hNode;
class hGraph {
    
private:
    int NUM_NODES;
    MatrixXi _adjMatrix;
    Eigen::VectorXi _degVector;
    int _eulerChar;
    double _dimension;
    double _hamiltonian;
    std::vector <int> _numCliques;
    bool cliquesFound = false;
    
public:
    hGraph(int size, MatrixXi adjMatrix);
    hGraph(int size);
    hGraph();
    ~hGraph();
    void print(); //depricated
    void toStream(std::ostream &os) const;
    void toFile(std::ofstream &fs) const;
    void setMatrix(MatrixXi data);
    void setMatrix(int size, MatrixXi data);
    void setHamiltonian(double val);
    hGraph unitSphere(int node);
    double calcDimension();
    void calcEulerChar();

    void numCliques();


    double getDimension();
    int getDegree(int node);
    int getSize();
    int getEulerChar();
    bool isConnected(int row, int column);
    void countCliques();
    void cliqueSearch(std::vector<int> R, std::vector<int> P);
    double getHam();
    void accept(absHamiltonian &ham);
    friend std::ostream &operator << (std::ostream &os, const hGraph &rhs);
    friend std::ofstream &operator <<(std::ofstream &fs, const hGraph &rhs);
};

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

hGraph * readGraphFile(int &num);

void removeColumn(Eigen::MatrixXi& matrix, unsigned int colToRemove);
void removeRow(Eigen::MatrixXi& matrix, unsigned int rowToRemove);
hGraph randomGraph(int size);
unsigned long int factorial(unsigned int n);
unsigned long int binom(unsigned int n, unsigned int k);
hGraph kGraph(int size);





#endif /* _hGraph_h */
