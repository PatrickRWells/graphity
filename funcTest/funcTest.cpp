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
#include <fstream>
#include <chrono>




using namespace std;

bool getTF();



int main() {

    int size = 50;
    hGraph test(size);
    test = randomGraph(size, 0);
    test.setThreads(4);
    auto start = std::chrono::high_resolution_clock::now();
    start = std::chrono::high_resolution_clock::now();
    test.calcDimension();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::high_resolution_clock::now()-start);
    
    std::cout << test.getDimension() << std::endl;



}

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
