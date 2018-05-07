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
    std::ofstream outfile;
    outfile.open("test.csv");
    hGraph test = kGraph(50);
    outfile << test;

}

