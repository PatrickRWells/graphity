



#include <iostream>
#include <string>
#include "hamiltonians.h"
#include <cmath>

int main() {
    double TINV;
    int SIZE;
    
    std::cout << "How large of a graph would you like to use? ";
    std::cin >> SIZE;
    std::cout << std::endl << "Input the inverse temperature: ";
    std::cin >> TINV;
    
    std::ofstream file;
    file.open("output.csv");
    
    
    hGraph * graph = new hGraph(SIZE);
    (*graph) = randomGraph(SIZE);

    std::function<void(hGraph&)> simFunction;
    std::function<double(hGraph&,int, int)> simPartial;
    
	simFunction = basicSquareHam;
	simPartial = basicSquarePartial;
    
    std::random_device rd; //c++11 random number generator
    std::mt19937 gen(rd());
    std::uniform_int_distribution<int> randNode(0, SIZE - 1);
    std::uniform_real_distribution<> probDis(0.0, 1.0);
    
    bool run = true;
    int nodeA = 0;
    int nodeB = 0;
    int sweeps = 0;
    int n = 0;
    int increases = 0;
    int swaps = 0;
    
    while(run) {
        
        while(nodeA == nodeB) {
            nodeA = randNode(gen);
            nodeB = randNode(gen);

        }

    
        simFunction(*graph);
        double hamInitial = graph->getHam();
        double hamDiff = simPartial(*graph, nodeA, nodeB);
        if(hamDiff <= 0) {
            graph->flipEdge(nodeA,nodeB);
            graph->acceptPartial(hamDiff);
            swaps++;
        }
        else {
            
            double prob = exp(TINV*(-hamDiff));
            if(probDis(gen) <= prob) {
                graph->flipEdge(nodeA,nodeB);
                graph->acceptPartial(hamDiff);
                increases++;
                swaps++;
            }

            
        }
        nodeA = 0;
        nodeB = 0;
        n++;
        if(n == 150) {
            n = 0;
            sweeps++;
            if(sweeps %100 == 0) {
                std::cout << sweeps << " sweeps performed" << std::endl;
            }

        }
        
        if(sweeps == 1000) {
            run = false;
        }
        
        

    }
    std::cout << "Graph outputed to file" << std::endl;
    file << *graph;
    std::cout << "Moves that increased energy: " << increases << std::endl;
    std::cout << "total swaps accepted: " << swaps <<std::endl;
    file.close();
    
    
}



