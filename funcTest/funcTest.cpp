//
//  funcTest.cpp
//  
//
//  Allows the user to test out functions.
//
//

#include <iostream>
#include "hGraph/hGraph.h"
#include "hamiltonians.h"
#include <algorithm>
#include <string>
using namespace std;


int main() {
    hGraph test = kGraph(20);
    basicSquareHam(test);
    double testA = partialBasicSquare(test, 1, 6);
    std::cout << test.getHam() << ", " << testA << std::endl;
    test.acceptPartial(testA);
    std::cout << test.getHam() << std::endl;
}


