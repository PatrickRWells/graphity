//
//  eulerChar.h
//
//
//  Created by Patrick on 11/17/17.
//
//  To use this eulerChar, start by changing all instances of "Template" to the name of your file (with correct capitalization)

#ifndef eulerChar_h
#define eulerChar_h

#include "absHamiltonian.h"
#include <cmath>

void EulerCharHam(hGraph &host, double source);
class EulerChar : public absHamiltonian {
    
private:
    double _result; //Result of the computation.
    int _partial;
    std::vector<int> nodeA;
    std::vector<int> nodeB;
    bool isPartial = false;
public:
    EulerChar(double source);
    EulerChar(std::vector<int> nodeA, std::vector<int> nodeB);
    double result();
    void calculate(hGraph &host); //Don't mess with these
    double getDifference() {
        return _partial;
    }

    
};

EulerChar::EulerChar(double source) : _result(0.0) { //unlikely that it should be changed, unless you have some background energy level
    sourceT = source;
    
}

EulerChar::EulerChar(std::vector<int> node1, std::vector<int> node2) : _result(0.0) { //unlikely that it should be changed, unless you have some background energy level
    nodeA = node1;
    nodeB = node2;
    isPartial = true;
    
}

double EulerChar::result() { //There should be no reason to edit this function;
    return _result;
}

void EulerChar::calculate(hGraph &host) { //This is where all the main calculation takes place.
    if(!isPartial) {
       _result = host.getEulerChar();   //
    }
    else {
        hGraph temp(host.getSize());
        temp = host;
        temp.flipEdge(nodeA,nodeB);
        EulerCharHam(temp, sourceT);
        double diff = temp.getHam() - host.getHam();
        _partial = diff;
    }

    
}

void EulerCharHam(hGraph &host, double source) { //do not edit this function except to change instances of the word "EulerChar"
                                 //For example, if your hamiltonian is named "basic" this function should be titled "basicHam"
    EulerChar Ham(source);
    host.accept(Ham);
    host.setHamiltonian(Ham.result());
    
}

double EulerCharPartial(hGraph &host, std::vector<int> nodeA, std::vector<int> nodeB) {
    EulerChar Ham(nodeA, nodeB);
    host.accept(Ham);
    return Ham.getDifference();
    
    
}


#endif /* basic2_h */
