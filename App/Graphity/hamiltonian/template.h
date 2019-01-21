//
//  template.h
//
//
//  Created by Patrick on 11/17/17.
//  While this template may be useful, it is no longer necessary to manually edit it to create new hamiltonians
//  Using the python script to do so is encouraged
//  To use this template, start by changing all instances of "Template" to the name of your file (with correct capitalization)

#ifndef template_h
#define template_h

#include "absHamiltonian.h"
#include <cmath>

void TemplateHam(hGraph &host);
class Template : public absHamiltonian {
    
private:
    double _result; //Use this attribute to store the result of full energy calculations.
    double _partial;//Use this attribute to store the result of partial energy calculations.
    std::vector<int> nodeA;
    std::vector<int> nodeB;
    int numEdges = 0;
    double sourceT = 0;
    bool isPartial = false;
public:
    Template();
    Template(std::vector<int> node1, std::vector<int> node2);
    double result();
    void calculate(hGraph &host); //Don't mess with these
    double getDifference() {
        return _partial;
    }
    
};

Template::Template() : _result(0.0) { //unlikely that it should be changed, unless you have some background energy level
    sourceT = SOURCE;
}

Template::Template(std::vector<int> node1, std::vector<int> node2) : _result(0.0) { //This constructor is used by the utility function when doing only a partial calculation.
    nodeA = node1;
    nodeB = node2;
    numEdges = node1.size();
    isPartial = true;
    
}


double Template::result() { //There should be no reason to edit this function;
    return _result;
}

void Template::calculate(hGraph &host) { //This is where all the main calculation takes place.
    int size = host.getSize();           //see hGraph documentation for available commands.
    
    if(isPartial) {                                //This function can take care of both full and partial calculations.
        hGraph temp(host.getSize()) = host;        //It is determined by the isPartial attribute, which is set  matically if the second constructor is used.
        for(int i = 0; i < numEdges; i++) {
            temp.flipEdge(nodeA[i], nodeB[i])'
        }
        _partial = 10;

    }
    else {
        for(int i = 0; i < size; i++) {
            _result += host.getDegree(i);   //
        }
    }
    
}

void TemplateHam(hGraph &host) { //do not edit this function except to change instances of the word "Template"
                                 //For example, if your hamiltonian is named "basic" this function should be titled "basicHam"
    Template Ham;
    host.accept(Ham);
    host.setHamiltonian(Ham.result());
    
}
//Note the difference between the partial hamiltonian and the full one. The full one sets the value in the graph object, the partial just returns a value
double TemplatePartial(hGraph &host, std::vector<int> nodeA, std::vector<int> nodeB) {    //Do not edit this function except to change instance of the word "Template"
    if(nodeA.size() != nodeB.size()) {
        std::cout << "Fatal error: vectors passed to a partial hamiltonian must contain the same number of elements" << std::endl;
        exit(2);
    }
    Template Ham(nodeA, nodeB);
    host.accept(Ham);
    return Ham.getDifference();
    
    
}

#endif /* basic2_h */
