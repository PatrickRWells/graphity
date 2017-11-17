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

#include "../lib/eigen/Eigen/Dense"

using Eigen::MatrixXi;
class absHamiltonian;

class hGraph {
    
private:
    int NUM_NODES;
    MatrixXi _adjMatrix;
    Eigen::VectorXi _degVector;
    double _eulerChar;
    double _hamiltonian;
    
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
    int getDegree(int node);
    int getSize();
    void accept(absHamiltonian &ham);
    friend std::ostream &operator << (std::ostream &os, const hGraph &rhs);
    friend std::ofstream &operator <<(std::ofstream &fs, const hGraph &rhs);
};

hGraph * readGraphFile(int &num);






#endif /* _hGraph_h */
