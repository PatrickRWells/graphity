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

using Eigen::MatrixXd;

class hGraph {
    
private:
    const int NUM_NODES;
    MatrixXd _adjMatrix;
    Eigen::VectorXd _degVector;
    double _eulerChar = 0;
    double _hamiltonian = 0;
    
public:
    hGraph(int size);
    ~hGraph();
    void print();
};




#endif /* _hGraph_h */
