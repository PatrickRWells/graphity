//
//  MonteCarlo.cpp
//  
//
//  Created by Patrick on 2/19/18.
//
//


#include <iostream>
#include <string>
#include "hamiltonians.h"
#include <cmath>

int main() {
    const double K_BOLTZ = 1.3806e-23;
    const int SIZE = 10;
    hGraph * graph = new hGraph(SIZE);
    (*graph) = randomGraph(SIZE);

    std::function<void(hGraph&)> simFunction;
	simFunction = basicSquareHam;
    
    std::random_device rd; //c++11 random number generator
    std::mt19937 gen(rd());
    std::uniform_int_distribution<int> randNode(0, SIZE - 1);
    std::uniform_real_distribution<> probDis(0.0, 1.0);
    
    
    hGraph * newGraph = new hGraph(SIZE);
    
    bool test = false;
    
    
    
    int nodeA = 0;
    int nodeB = 0;
    
    while(!test) {
        
        while(nodeA == nodeB) {
            nodeA = randNode(gen);
            nodeB = randNode(gen);

        }

    
        *newGraph = *graph;
        newGraph->flipEdge(nodeA, nodeB);

    
        simFunction(*graph);
        simFunction(*newGraph);
    
        double hamInitial = graph->getHam();
        double hamFinal = newGraph->getHam();
        
        /*if(hamFinal >= hamInitial) {
            continue;
        }*/
        
            double prob = (exp(hamFinal))/(exp(hamInitial));
            double probA = probDis(gen);
        std::cout << prob << "," << probA << std::endl;
            
                test = true;

    

    
    
       /* if(hamFinal < hamInitial) {
        }*/

    }
    
    std::cout << *graph << std::endl;
    std::cout <<  *newGraph << std::endl;

    

    
    
}



