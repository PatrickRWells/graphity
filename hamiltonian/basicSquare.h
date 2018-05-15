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
    std::vector<int> nodeA;
    std::vector<int> nodeB;
    int numEdges = 0;
    double source;
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

basicSquare::basicSquare() : _result(0.0) {
    
}

basicSquare::basicSquare(std::vector<int> node1, std::vector<int> node2) : _result(0.0) { //unlikely that it should be changed, unless you have some background energy level
    //std::cout << "Input value of source term: ";
    //std::cin >> source;
    source = -0.1;
    nodeA = node1;
    nodeB = node2;
    numEdges = node1.size();
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
                    _result += (diff*diff)*0.1;
                
            }
        }
        
        if(source != 0) {
            int sum = 0;
            for(int i = 0; i < size; i++) {
                sum += host.getDegree(i);
            
            }
            sum /= 2;
            sum *= source;
            _result += sum;
        }
    }
    else {
        hGraph temp(host.getSize());
        temp = host;
        for(int i = 0; i < numEdges; i++) {
            temp.flipEdge(nodeA[i],nodeB[i]);

        }
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

#endif /* basic2_h */
