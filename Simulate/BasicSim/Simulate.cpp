//
//  Simulate.cpp
//  
//
//  Created by Patrick on 11/2/17.
//
//


#include <iostream>
#include <string>
#include "hamiltonians.h"

int main() {

    hGraph * graphData;
    int numRead;
    std::function<void(hGraph&)> simFunction;
    
	simFunction = basicSquareHam;
    
    hList * first;
    hList * second;
    hList *t;
    first = new hList;
    second = new hList;
    graphData = readGraphFile(numRead);
    
    std::cout << "Calculating hamiltonians..." << std::endl;
    
    for(int i = 0; i < numRead; i++) {
        simFunction(graphData[i]);
        
        
        if((graphData[i].getHam() < first->getEnergy()) || first->isEmpty()) {
            t = second;
            second = first;
            first = t;
            first->clear();
            first->prepend(graphData + i);
            
        }
        
        else if(graphData[i].getHam() == first->getEnergy()) {
            first->prepend(graphData + i);

        }

        
        else if((graphData[i].getHam() < second->getEnergy()) || second->isEmpty()) {
            second->clear();
            second->prepend(graphData + i);

        }
        
        
        else if(graphData[i].getHam() == second->getEnergy()){
            second->prepend(graphData + i);

        }
        
        
    }
    std::cout << "Output file will include lowest energy graphs, follow by those in the first excited state." <<std::endl;
    std::cout << "Please enter an output file name: ";
    std::string outFile;
    std::cin >> outFile;
    outFile += ".csv";
    std::ofstream ofile;
    ofile.open(outFile);
    ofile << *first;
    ofile << *second;
    ofile.close();
    
    delete [] graphData;
    
    
}



