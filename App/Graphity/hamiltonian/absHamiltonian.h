//
//  absHamiltonian.h
//  
//
//  Created by Patrick on 10/31/17.
//
//


#ifndef absHamiltonian_h
#define absHamiltonian_h

class hGraph;

class absHamiltonian {
public:
    virtual void calculate(hGraph &host) = 0;    
};


#endif /* absHamiltonian_h */

