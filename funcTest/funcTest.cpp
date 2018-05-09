//
//  funcTest.cpp
//  
//
//  Allows the user to test out functions. The makefile links everything that other modules in this project might use, including Root and Eigen.
//
//

#include <iostream>
#include "hGraph/hGraph.h"
#include "hamiltonians.h"
#include <algorithm>
#include <string>
#include <chrono>
using namespace std;


int main() {
    
    hGraph test = compGraph(10);
    std::cout << test << std::endl;
    
}
