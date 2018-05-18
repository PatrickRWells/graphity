//
//  funcTest.cpp
//  
//
//  Allows the user to test out functions. The makefile links everything that other modules in this project might use, including Root and Eigen.
//
//

#include <iostream>
#include "hGraph/hGraph.h"
#include "graphics/graphingUtil.hpp"
#include "hamiltonians.h"
#include <algorithm>
#include <string>
#include <chrono>
#include <string>
#include <iomanip>

using namespace std;


int main() {
    int numRead = 0;
    hGraph * test = readGraphFile(numRead);
    for(int i = 0; i < numRead; i++) {
        std::cout << test[i] << std::endl;
        
    }
    

    
}
