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

using namespace std;


int main() {
    
    hGraph test = randomGraph(20);
    std::vector<double> dimen = test.getSpectralDimen();
    for(int i = 0; i < dimen.size(); i++ ) {
        std::cout << dimen[i] << std::endl;
    }
    std::cout << dimen.size() << std::endl;


    
}
