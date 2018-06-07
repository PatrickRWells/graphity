//
//  curveDiff.h
//
//
//  Created by Patrick on 11/17/17.
//
//  To use this curveDiff, start by changing all instances of "curveDiff" to the name of your file (with correct capitalization)

#ifndef curveDiff_h
#define curveDiff_h

#include "absHamiltonian.h"
#include <cmath>

void curveDiffHam(hGraph &host);
class curveDiff : public absHamiltonian {
    
private:
    double _result; //Use this attribute to store the result of full energy calculations.
    double _partial;//Use this attribute to store the result of partial energy calculations.
    std::vector<int> nodeA;
    std::vector<int> nodeB;
    int numEdges = 0;
    bool isPartial = false;
    double sourceT = 0;
public:
    curveDiff();
    curveDiff(std::vector<int> node1, std::vector<int> node2);
    double result();
    void calculate(hGraph &host); //Don't mess with these
    double getDifference() {
        return _partial;
    }
    
};

curveDiff::curveDiff() : _result(0.0) { //unlikely that it should be changed, unless you have some background energy level
    sourceT = SOURCE;
}

curveDiff::curveDiff(std::vector<int> node1, std::vector<int> node2) : _result(0.0) { //This constructor is used by the utility function when doing only a partial calculation.
    nodeA = node1;
    nodeB = node2;
    numEdges = node1.size();
    isPartial = true;
    
}


double curveDiff::result() { //There should be no reason to edit this function;
    return _result;
}

void curveDiff::calculate(hGraph &host) { //This is where all the main calculation takes place.
    if(!isPartial) {
    
        int size = host.getSize();
        double data[size];
        double sum = 0;
        for(int i = 0; i < size; i++) {
            double val = host.curvatureAt(i);
            data[i] = val;
            sum += val;
            
        }
        sum /= size;
        double sum2 = 0;
        for(int i = 0; i < size; i++) {
            
            sum2 += (data[i] - sum)*(data[i] - sum);
            
        }
        _result = sum2;
        if(sourceT != 0) {
            for(int i = 0; i < size; i++) {
                _result += sourceT*host.getDegree(i);
            }
            
        }
    }
    
    else {
        hGraph temp(host.getSize());
        temp = host;
        temp.flipEdge(nodeA,nodeB);
        curveDiffHam(temp);
        double diff = temp.getHam() - host.getHam();
        _partial = diff;
        
    }
    

    
}

void curveDiffHam(hGraph &host) { //do not edit this function except to change instances of the word "curveDiff"
                                 //For example, if your hamiltonian is named "basic" this function should be titled "basicHam"
    curveDiff Ham;
    host.accept(Ham);
    host.setHamiltonian(Ham.result());
    
}
//Note the difference between the partial hamiltonian and the full one. The full one sets the value in the graph object, the partial just returns a value
double curveDiffPartial(hGraph &host, std::vector<int> nodeA, std::vector<int> nodeB) {    //Do not edit this function except to change instance of the word "curveDiff"
    if(nodeA.size() != nodeB.size()) {
        std::cout << "Fatal error: vectors passed to a partial hamiltonian must contain the same number of elements" << std::endl;
        exit(2);
    }
    curveDiff Ham(nodeA, nodeB);
    host.accept(Ham);
    return Ham.getDifference();
    
    
}

#endif /* basic2_h */
