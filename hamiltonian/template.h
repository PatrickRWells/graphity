//
//  template.h
//
//
//  Created by Patrick on 11/17/17.
//
//  To use this template, start by changing all instances of "Template" to the name of your file (with correct capitalization)

#ifndef template_h
#define template_h

#include "absHamiltonian.h"
#include <cmath>

void TemplateHam(hGraph &host);
class Template : public absHamiltonian {
    
private:
    double _result; //Result of the computation.
public:
    Template();
    double result();
    void calculate(hGraph &host); //Don't mess with these
    
    
};

Template::Template() : _result(0.0) { //unlikely that it should be changed, unless you have some background energy level
    
}

double Template::result() { //There should be no reason to edit this function;
    return _result;
}

void Template::calculate(hGraph &host) { //This is where all the main calculation takes place.
    int size = host.getSize();           //see hGraph documentation for available commands.
    for(int i = 0; i < size; i++) {
        _result += host.getDegree(i);   //
    }
    
}

void TemplateHam(hGraph &host) { //do not edit this function except to change instances of the word "Template"
                                 //For example, if your hamiltonian is named "basic" this function should be titled "basicHam"
    Template Ham;
    host.accept(Ham);
    host.setHamiltonian(Ham.result());
    
}

#endif /* basic2_h */
