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
public:
    basicSquare();
    double result();
    void calculate(hGraph &host); //Don't mess with these
    
    
};

basicSquare::basicSquare() : _result(0.0) { //unlikely that it should be changed, unless you have some background energy level
    
}

double basicSquare::result() { //just leave the function iteslf as-is
    return _result;
}

void basicSquare::calculate(hGraph &host) { //This is where all the main calculation takes place
    int size = host.getSize();
    for(int i = 0; i < size; i++) {
        for(int j = i+1; j < size; j++) {
            double diff = host.getDegree(i) - host.getDegree(j);
            _result += diff*diff;
            
        }
    }
    
}

void basicSquareHam(hGraph &host) {
    basicSquare Ham;
    host.accept(Ham);
    host.setHamiltonian(Ham.result());
    
}

#endif /* basic2_h */
