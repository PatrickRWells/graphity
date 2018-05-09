//
//  basicSquare.h
//
//
//  Created by Patrick on 11/17/17.
//
//  To use this basicSquare, start by changing all instances of "basicSquare" to the name of your file

#ifndef basicSquare_h
#define basicSquare_h

#include "absHamiltonian.h"
#include <cmath>

void basicSquareHam(hGraph &host);
class basicSquare : public absHamiltonian {
    
private:
    int _result; //Result of the computation. Should not be changed.
    int _partial;
    int nodeA;
    int nodeB;
    bool isPartial = false;
public:
    basicSquare();
    basicSquare(int node1, int node2);
    double result();
    void calculate(hGraph &host); //Don't mess with these
    double getDifference() {
        return _partial;
    }
    
};

basicSquare::basicSquare() : _result(0.0) {
    
}

basicSquare::basicSquare(int node1, int node2) : _result(0.0) { //unlikely that it should be changed, unless you have some background energy level
    nodeA = node1;
    nodeB = node2;
    isPartial = true;
    
}

double basicSquare::result() { //just leave the function iteslf as-is
    return _result;
}

void basicSquare::calculate(hGraph &host) { //This is where all the main calculation takes place
    if(!isPartial) {
        int size = host.getSize();
        for(int i = 0; i < size; i++) {
            for(int j = i+1; j < size; j++) {
                    double diff = host.getDegree(i) - host.getDegree(j);
                    _result += diff*diff;
                
            }
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

double basicSquarePartial(hGraph &host, int nodeA, int nodeB) {
    basicSquare Ham(nodeA, nodeB);
    host.accept(Ham);
    return Ham.getDifference();
    
    
}

#endif /* basic2_h */
