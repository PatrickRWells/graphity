//
//  expectationVals.cpp
//  
//
//  Created by Patrick on 5/21/18.
//

#include <stdio.h>
#include "hGraph/hGraph.h"
#include "hamiltonians.h"

int main() {
    
    
    int numRead;
    hGraph ** graphs;
    readGraphFile(&graphs, numRead);
    std::cout << numRead << " read." << std::endl;
    std::cout << *graphs[numRead -1] << std::endl;

    

}
