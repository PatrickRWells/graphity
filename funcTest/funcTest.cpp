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
#include <chrono>
using namespace std;


int main() {
    
    hGraph test = compGraph(10);
    std::cout << test << std::endl;
    
}


/*
 
 
 Calculating shortest cycle at node k
 
 Find all neighbors of k and keep track of all walks. Obviously none are cycles
 Make another step from each of the resultant paths
 
 Check if any of the last items in the three-walks is equal to the last node in the two-walk. If so, the shortest cycle is 3.
 
 Then go to four-walks. If any of the last items in the 4-walks is the last item in one of the two-walks, the shortest cycle is four.
 If not, check all four walks against eachother. If two cycles have the same end points, and are not identical, then the shortest cycle is four.
 
 At each iteration, must check all previous walks and see if two walks terminate on the same node. However if they contain a shared internal node, then it is not a loop.  
 
 
 
 
 
 */
