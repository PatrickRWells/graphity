//
//  basicSquare.h
//
//
//  Created by Patrick on 11/17/17.
//
//  A simple hamiltonian that can be used to examine the properties of the model

#ifndef basicSquare_h
#define basicSquare_h

#include "absHamiltonian.h"
#include <cmath>


void basicSquareHam(hGraph &host, double source);
class basicSquare : public absHamiltonian {
    
private:
    double _result; //Result of the computation. Should not be changed.
    double _partial;
    std::vector<int> nodeA;
    std::vector<int> nodeB;
    int numEdges = 0;
    bool isPartial = false;
public:
    basicSquare(double source);
    basicSquare(std::vector<int> node1, std::vector<int> node2);
    double result();
    void calculate(hGraph &host); //Don't mess with these
    double getDifference() {
        return _partial;
    }
    
};


basicSquare::basicSquare(double source) : _result(0.0) {
    sourceT = source;
}

basicSquare::basicSquare(std::vector<int> node1, std::vector<int> node2) : _result(0.0) { //unlikely that it should be changed, unless you have some background energy level
    //std::cout << "Input value of source term: ";
    //std::cin >> source;
    nodeA = node1;
    nodeB = node2;
    numEdges = node1.size();
    isPartial = true;
    
}

double basicSquare::result() { //just leave the function iteslf as-is
    return _result;
}

void basicSquare::calculate(hGraph &host) { //This is where all the main calculation takes place
    
    if(!isPartial) { //Full calculation
        int size = host.getSize();
        for(int i = 0; i < size; i++) {
            for(int j = i+1; j < size; j++) {
                double diff = host.getDegree(i) - host.getDegree(j);
                _result += (diff*diff); //The sum of the square of the difference between the degree of every pair of nodes on the graphs
            
            }
        }
        _result /= (2*host.getSize()); //Prefactor of 1/2N
        
        if(sourceT != 0) {
            double sum = 0;
            for(int i = 0; i < size; i++) {
                sum += host.getDegree(i); //A source term that tunes the prefered value of the node degree
            
            }
            sum *= sourceT;
            _result += sum;
        }
    }
    
    
    
    else { //A partial calculation just creates a new graph with the nodes flipped and finds the value of the hamiltonian.
        hGraph temp(host.getSize());
        temp = host;
        temp.flipEdge(nodeA,nodeB);

        basicSquareHam(temp, sourceT);
        double diff = temp.getHam() - host.getHam();
        _partial = diff;
    }
    
}

void basicSquareHam(hGraph &host, double source) {
    basicSquare Ham(source);
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
