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
    double sourceT = 0;
    bool isPartial = false;
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
        host.getEulerChar(); //Has to calculate the curvature at every point anyway and is multithreaded.
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

void curveDiffHam(hGraph &host) {
                                
    curveDiff Ham;
    host.accept(Ham);
    host.setHamiltonian(Ham.result());
    
}

double curveDiffPartial(hGraph &host, std::vector<int> nodeA, std::vector<int> nodeB) { 
    if(nodeA.size() != nodeB.size()) {
        std::cout << "Fatal error: vectors passed to a partial hamiltonian must contain the same number of elements" << std::endl;
        exit(2);
    }
    curveDiff Ham(nodeA, nodeB);
    host.accept(Ham);
    return Ham.getDifference();
    
    
}

#endif /* curveDiff_h */
