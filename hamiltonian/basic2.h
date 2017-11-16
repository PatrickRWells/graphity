//
//  basic.h
//
//
//  Created by Patrick on 10/31/17.
//
//

#ifndef basic_h
#define basic_h

#include "absHamiltonian.h"
#include <cmath>

class basic : public absHamiltonian {
    
private:
    int _result;
public:
    basic();
    double result();
    void calculate(hGraph &host);
    
    
};

basic::basic() : _result(0.0) {
    
}

double basic::result() {
    return _result;
}

void basic::calculate(hGraph &host) {
    int size = host.getSize();
    std::cout << size << std::endl;
    for(int i = 0; i < size; i++) {
        _result += host.getDegree(i);
    }
    std::cout << _result << std::endl;
    
}


double basicH(hGraph host) {
    basic basicCalc;
    host.accept(basicCalc);
    return basicCalc.result();
    
}


#endif /* basic_h */
