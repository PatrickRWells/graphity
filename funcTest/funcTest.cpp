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
#include <chrono>

#define TIMING

#ifdef TIMING
#define INIT_TIMER auto start = std::chrono::high_resolution_clock::now();
#define START_TIMER  start = std::chrono::high_resolution_clock::now();
#define STOP_TIMER(name)  std::cout << "RUNTIME of " << name << ": " << \
std::chrono::duration_cast<std::chrono::milliseconds>( \
std::chrono::high_resolution_clock::now()-start \
).count() << " ms " << std::endl;
#else
#define INIT_TIMER
#define START_TIMER
#define STOP_TIMER(name)
#endif


using namespace std;

bool getTF() {
    while(true) {
        char c = getchar();
        if (toupper(c) == 'Y') {
            std::cin.clear();
            std::cin.ignore(100, '\n');
            return true;
            
        }
        else if(toupper(c) == 'N') {
            std::cin.clear();
            std::cin.ignore(100, '\n');
            return false;
            
        }
        else {
            std::cout << "Invalid input." << std::endl;
            std::cout << "Please enter a 'y' or an 'n': ";
            std::cin.clear();
            std::cin.ignore(100, '\n');
            
        }
    }
}

int main() {
    INIT_TIMER;
    int size = 40;
    hGraph test(size);
    test = randomGraph(size);
    test.setThreads(4);
    
    std::cout << "Old dimension calculation: ";
    START_TIMER
    test.oldCalcDimension();
    STOP_TIMER("Old dimension calculation")
    
    std::cout << "New dimension calculation: ";
    START_TIMER
    test.calcDimension();
    STOP_TIMER("New dimension calculation")

    std::cout << "Woud you like to output the graph? (y/n) ";
    if(getTF()) {
        cout << test;
    }
    //test.flipEdge(xVals, yVals);
    
}
