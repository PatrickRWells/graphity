//
//  basic.h
//
//
//  Created by Patrick on 10/31/17.
//
//

#ifndef basic2_h
#define basic2_h

#include "absHamiltonian.h"
#include <cmath>
void basic2Ham(hGraph &host);
class basic2 : public absHamiltonian {
    
private:
    int _result;
public:
    basic2();
    double result();
    void calculate(hGraph &host);
    
    
};

basic2::basic2() : _result(0.0) {
    
}

double basic2::result() {
    return _result;
}

void basic2::calculate(hGraph &host) {
    std::cout << "test" << std::endl;
    int size = host.getSize();
    std::cout << size << std::endl;
    for(int i = 0; i < size; i++) {
        _result += host.getDegree(i);
    }
    std::cout << _result << std::endl;
    
}

void basic2Ham(hGraph &host) {
    basic2 basic2Ham;
    host.accept(basic2Ham);
    std::cout << basic2Ham.result() << std::endl;
    host.setHamiltonian(basic2Ham.result());
    
}

#endif /* basic2_h */
