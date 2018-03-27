//
//  stats.cpp
//  
//
//  Created by Patrick on 3/13/18.
//

#include <iostream>
#include "hGraph/hGraph.h"
#include "hamiltonians.h"
#include <cmath>
#include <iomanip>
#include <thread>
#include <atomic>

const int SIZE = 100;
const int NUMGRAPHS = 1000;
const int NUMTESTS = 1000;
std::atomic<int> completed(0);

void calculate(double **data, int numgraphs, int numtests, int size, int offset, std::function<void(hGraph&)> simulate);


int main() {
    std::function<void(hGraph&)> simFunction;
    simFunction = basicSquareHam;

    std::cout << "Simulating..." << std::endl;
    double ** data;
    data = new double*[NUMGRAPHS];
    for(int i = 0; i < NUMGRAPHS; i++) {
        data[i] = new double[NUMTESTS+1];

    }
    std::thread t1(calculate, data, 250, NUMTESTS, SIZE, 0, simFunction);
    std::thread t2(calculate, data, 250, NUMTESTS, SIZE, 250, simFunction);
    std::thread t3(calculate, data, 250, NUMTESTS, SIZE, 500, simFunction);
    calculate(data, 250, NUMTESTS, SIZE, 750, simFunction);
    
    t1.join();
    t2.join();
    t3.join();
    double sum = 0;
    for(int i = 0; i < NUMGRAPHS; i++) {
        for(int j = 0; j < NUMTESTS; j++) {
            sum += data[i][j];
        }
    }
    double avg = sum/(NUMGRAPHS*NUMTESTS);
    sum = 0;
    for(int i = 0; i < NUMGRAPHS; i++) {
        for(int j = 0; j < NUMTESTS; j++) {
            double val = data[i][j] - avg;
            sum += val*val;
        }
    }
    
    sum = sum/(NUMGRAPHS*NUMTESTS);
    double stddev = sqrt(sum);

    double sum3;
    for(int i = 0; i < NUMGRAPHS; i++) {
        sum3 += data[i][NUMTESTS];
    }
    double avg3 = sum3/NUMGRAPHS;
    sum3 = 0;
    for(int i = 0; i < NUMGRAPHS; i++) {
            double val = data[i][NUMTESTS] - avg3;
            sum3 += val*val;

    }

    double stddev3 = sqrt(sum3/NUMGRAPHS);
    std::cout << "Complete               " << std::endl;
    
    std:: cout << "Average: " << avg << std::endl;
    std::cout << "Standard Deviation: " << stddev << std::endl;
    std::cout << "Average energy: " << avg3 << std::endl;
    std::cout << "Standard Deviation: " << stddev3 << std::endl;

}

void calculate(double ** data, int numgraphs, int numtests, int size, int offset, std::function<void(hGraph&)> simulate) {
    
    std::random_device rd; //c++11 random number generator
    std::mt19937 gen(rd());
    std::uniform_int_distribution<int> randNode(0, size - 1);
    int nodeA = 0;
    int nodeB = 0;
    
    for(int i = offset; i < numgraphs + offset; i++) {
        hGraph * graph = new hGraph(size);
        (*graph) = randomGraph(size);
        simulate(*graph);
        data[i][NUMTESTS] = graph->getHam();
        for (int j = 0; j < numtests; j++) {
            data[i][j] = -(graph->getHam());
            hGraph * newGraph = new hGraph(size);
            *newGraph = *graph;
            do {
                nodeA = randNode(gen);
                nodeB = randNode(gen);
                
            } while (nodeA == nodeB);
            newGraph->flipEdge(nodeA, nodeB);
            simulate(*newGraph);
            data[i][j] += newGraph->getHam();
            delete newGraph;
        }
        delete graph;
        completed++;
        if(completed.load() % 25 == 0) {
            std::cout << completed.load() << " graphs tested" << std::endl;
        }

        //std:: cout << " Progress - " << std::setprecision(3) <<  100*i/float(NUMGRAPHS) << " %\r         ";
        //std::cout.flush();
        
    }
    
    
    
}
