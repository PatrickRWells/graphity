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
    hGraph(int size);
    ~hGraph();
    void print();
    void setMatrix(MatrixXi data);
    void setVector(Eigen::VectorXi data);
};




#endif /* _hGraph_h */
