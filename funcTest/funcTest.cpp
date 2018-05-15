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
    test.flipEdge(0,1);
    test.flipEdge(2,3);
    test.flipEdge(4,5);
    test.flipEdge(6,7);
    test.flipEdge(8,9);
    test.flipEdge(10,11);
    test.flipEdge(12,13);
    std::cout << test.getDimension() << std::endl;


    
}
