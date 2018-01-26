//
//  funcTest.cpp
//  
//
//  Allows the user to test out functions.
//
//

#include <iostream>
#include "hGraph/hGraph.h"
#include <algorithm>
#include <string>
using namespace std;


int main() {
    hGraph test = kGraph(20);
    test.calcEulerChar();
    std::cout << test ;

    
    
}


