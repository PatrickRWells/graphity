#ifndef basicSquare_h
#define basicSquare_h

#include "absHamiltonian.h"
#include <cmath>

void basicSquareHam(hGraph &host);
class basicSquare : public absHamiltonian {
    
private:
    double _result; //Use this attribute to store the result of full energy calculations.
    double _partial;//Use this attribute to store the result of partial energy calculations.
    std::vector<int> nodeA;
    std::vector<int> nodeB;
    int numEdges = 0;
    double sourceT = 0;
    bool isPartial = false;
public:
    basicSquare();
    basicSquare(std::vector<int> node1, std::vector<int> node2);
    double result();
    void calculate(hGraph &host); //Don't mess with these
    double getDifference() {
        return _partial;
    }
    
};

basicSquare::basicSquare() : _result(0.0) { //unlikely that it should be changed, unless you have some background energy level
    sourceT = SOURCE;
}

basicSquare::basicSquare(std::vector<int> node1, std::vector<int> node2) : _result(0.0) { //This constructor is used by the utility function when doing only a partial calculation.
    nodeA = node1;
    nodeB = node2;
    numEdges = node1.size();
    isPartial = true;
    
}


double basicSquare::result() { //There should be no reason to edit this function;
    return _result;
}

void basicSquare::calculate(hGraph &host) { //This is where all the main calculation takes place.    
    if(!isPartial) {
        int size = host.getSize();
        for(int i = 0; i < size; i++) {
            for(int j = i+1; j < size; j++) {
                double diff = host.getDegree(i) - host.getDegree(j);
                _result += (diff*diff);
            
            }
        }
        _result /= (2*host.getSize());
        
        if(sourceT != 0) {
            double sum = 0;
            for(int i = 0; i < size; i++) {
                sum += host.getDegree(i);
            
            }
            sum *= sourceT;
            _result += sum;
        }
    }
    
    
    
    else {
        hGraph temp(host.getSize());
        temp = host;
        temp.flipEdge(nodeA,nodeB);

        basicSquareHam(temp);
        double diff = temp.getHam() - host.getHam();
        _partial = diff;
    }
}

void basicSquareHam(hGraph &host) {
                                
    basicSquare Ham;
    host.accept(Ham);
    host.setHamiltonian(Ham.result());
    
}

double basicSquarePartial(hGraph &host, std::vector<int> nodeA, std::vector<int> nodeB) { 
    if(nodeA.size() != nodeB.size()) {
        std::cout << "Fatal error: vectors passed to a partial hamiltonian must contain the same number of elements" << std::endl;
        exit(2);
    }
    basicSquare Ham(nodeA, nodeB);
    host.accept(Ham);
    return Ham.getDifference();
    
    
}

#endif /* basicSquare_h */
