//
//  Simulate.cpp
//  
//
//  Created by Patrick on 11/2/17.
//
//


#include <iostream>
#include <string>
#include "hamiltonians.h"

int main() {

    hGraph * graphData;
    int numRead;
    std::function<void(hGraph&)> simFunction;
    
	simFunction = basicHam;
    
    
    graphData = readGraphFile(numRead);
    std::cout << graphData[0] << std::endl;
    simFunction(graphData[0]);
    std::cout << graphData[0] << std::endl;

    
    
}



