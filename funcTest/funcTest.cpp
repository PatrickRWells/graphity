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
    
    std::vector<std::vector<int>> xVals;
    std::vector<std::vector<int>> yVals;
    
    std::vector<int> tempX;
    std::vector<int> tempY;
    for(int i = 0; i < 10; i++) {
        tempX.push_back(i);
        tempY.push_back(2*i);
        
    }
    
    xVals.push_back(tempX);
    yVals.push_back(tempY);
    
    tempX.clear();
    tempY.clear();
    
    for(int i = 0; i < 10; i++) {
        tempX.push_back(i);
        tempY.push_back(-2*i + 20);
        
    }
    
    xVals.push_back(tempX);
    yVals.push_back(tempY);
    
    tempX.clear();
    tempY.clear();
    
    for(int i = 0; i < 10; i++) {
        tempX.push_back(i);
        tempY.push_back(10);
        
    }
    
    xVals.push_back(tempX);
    yVals.push_back(tempY);
    
    drawMultiGraph(xVals, yVals);


    


    
}
