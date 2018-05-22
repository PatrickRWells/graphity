//
//  expectationVals.cpp
//  
//
//  Created by Patrick on 5/21/18.
//

double TINV = 1;
#include <stdio.h>
#include <cmath>
#include "hGraph/hGraph.h"
#include "hamiltonians.h"

int main() {
    
    SOURCE = -10;
    
    
    int numRead;
    hGraph ** graphs;
    readGraphFile(&graphs, numRead);
    double partition = 0.0;
    for(int i = 0; i < numRead; i++) {
        basicSquareHam(*graphs[i]);
        partition += exp(-TINV*graphs[i]->getHam());
    }
    double sum = 0.0;
    for(int i = 0; i < numRead; i++) {
        sum+= (graphs[i]->avgDegree())*(exp(-TINV*graphs[i]->getHam()))/partition;
        
        
    }
    
    std::cout << "Average degree expecation value: " << sum << std::endl;
    delete [] graphs;
    
    
    
    

    

}
