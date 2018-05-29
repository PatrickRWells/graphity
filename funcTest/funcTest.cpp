//
//  funcTest.cpp
//  
//
//  Allows the user to test out functions. The makefile links everything that other modules in this project might use, including Root and Eigen.
//
//

#include <iostream>
#include "hGraph/hGraph.h"
#include "graphics/graphUtil/graphingUtil.hpp"
#include "hamiltonians.h"
#include <algorithm>
#include <string>
#include <chrono>
#include <string>
#include <iomanip>

using namespace std;


int main() {
    int size = 8;
    hGraph test(size);
    test = compGraph(size);
    test.flipEdge(0,1);
    test.flipEdge(0,2);
    test.flipEdge(0,3);

    test.flipEdge(4,5);
    test.flipEdge(4,6);
    test.flipEdge(1,5);
    std::cout << test.getDimension() << std::endl;

    

    //test.flipEdge(xVals, yVals);
    
}
