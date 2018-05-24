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
    hGraph test(10);
    test = randomGraph(10);
    SOURCE = 0;

    //test.flipEdge(xVals, yVals);
    basicSquareHam(test);
    std::cout << test << std::endl;
    
}
