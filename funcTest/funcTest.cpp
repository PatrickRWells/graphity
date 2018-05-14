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
    
    hGraph test = compGraph(20);
    test.flipEdge(1,2);
    test.setThreads(4);
    std::cout << test.getDimension() << std::endl;


    
}
