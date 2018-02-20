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
    int _eulerChar = 0;
    double _dimension;
    double _hamiltonian;
    bool cliquesFound = false;


    MatrixXi _adjMatrix;
    Eigen::VectorXi _degVector;
    std::vector <int> _numCliques;
    
    void toStream(std::ostream &os) const;
    void toFile(std::ofstream &fs) const;
    void accept(absHamiltonian &ham);


    
public:
    hGraph(int size, MatrixXi adjMatrix);
    hGraph(int size);
    hGraph();
    ~hGraph();
    void print(); //depricated
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

hGraph * readGraphFile(int &num);

void removeColumn(Eigen::MatrixXi& matrix, unsigned int colToRemove);
void removeRow(Eigen::MatrixXi& matrix, unsigned int rowToRemove);
hGraph randomGraph(int size);
hGraph kGraph(int size);





#endif /* _hGraph_h */
