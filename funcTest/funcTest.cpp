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
    
    hGraph graph(size);
    graph = zeroGraph(size);
    //graph = compGraph(size);

    std::cout << graph;

    while(true) {
        std::vector<int> initialDimen = graph.fractionalDimen();
        char nodeA;
        char nodeB;
        std::cin >> nodeA;
        std::cin >> nodeB;
        graph.flipEdge(nodeA - '0', nodeB - '0');
        std::vector<int> finalDimen = graph.fractionalDimen();
        fractReduce(finalDimen);

        initialDimen[0] *= -1;
        std::cout << graph;
        std::cout << "Fractional dimensionality: " << finalDimen[0] << "/" << finalDimen[1] << std::endl;
        initialDimen = fractionAdd(initialDimen, finalDimen);
        std::cout << "Fractional change in dimensionality: " << initialDimen[0] << "/" << initialDimen[1] << std::endl;
        
    }
    std::vector<int> temp = graph.fractionalDimen();
    fractReduce(temp);
    std::cout << temp[0] << "/" << temp[1] << std::endl;

    //test.flipEdge(xVals, yVals);
    
}
