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


class hGraph {
    
private:
    const int NUM_NODES;
    MatrixXi _adjMatrix;
    Eigen::VectorXi _degVector;
    double _eulerChar = 0;
    double _hamiltonian = 0;
    
public:
    hGraph(int size, MatrixXi adjMatrix);
    hGraph(int size);
    ~hGraph();
    void print();
    void toStream(std::ostream &os) const;
    void toFile(std::ofstream &fs) const;
    void setMatrix(MatrixXi data);
    friend std::ostream &operator << (std::ostream &os, const hGraph &rhs);
    friend std::ofstream &operator <<(std::ofstream &fs, const hGraph &rhs);
};






#endif /* _hGraph_h */
