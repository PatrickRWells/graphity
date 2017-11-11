//
//  Simulate.cpp
//  
//
//  Created by Patrick on 11/2/17.
//
//

#include <iostream>
#include "hGraph/hGraph.h"
#include "hamiltonians.h"
#include <cmath>

int main() {
    
    char cSize;
    std::ifstream input;
    std::ofstream output;
    input.open("../test.csv");
    input >> cSize;
    int size = cSize - '0';
    int num = size*(size-1)/2;
    num = pow(2, num) - 1;
    
    }
    
    
    
    
    
}


