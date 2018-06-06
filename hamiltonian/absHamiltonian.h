//
//  absHamiltonian.h
//  
//
//  Created by Patrick on 10/31/17.
//
//


#ifndef absHamiltonian_h
#define absHamiltonian_h
#include "hGraph/hGraph.h"

class absHamiltonian {
protected:
    static double sourceT;
public:
    virtual void calculate(hGraph &host) = 0;
};

double absHamiltonian::sourceT = 0;

#endif /* absHamiltonian_h */

