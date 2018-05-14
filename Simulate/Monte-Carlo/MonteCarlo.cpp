#include <iostream>
#include <string>
#include "hamiltonians.h"
#include "graphics/graphingUtil.hpp"
#include <cmath>
#include <vector>
#include <thread>

void monteCarlo(hGraph * graph, std::vector<double> * energy, std::vector<double> * dimensions, bool progress, std::string descriptor);
void correlationFn(std::vector<double> * data, std::vector <double> * output); //Simulation itself must be a function in order to be paralelized.

int NUM_CORES = 4; //Number of cores in the CPU
double TINV;    //Beta
int SIZE;       //Number of nodes in the graphs
int maxSweeps;  //Number of sweeps that will be performed


//------Simulation Functions------//
std::function<void(hGraph&)> simFunction;
std::function<double(hGraph&,std::vector<int>,std::vector<int>)> simPartial;
//All these variables are declared to be global so that they can be used by multiple threads.



int main() {
    std::ofstream corrOut("correlation.csv");
    std::ofstream("dimensionality.csv");
    //Sets simulation functions
    simFunction = basicSquareHam;
    simPartial = basicSquarePartial;

    //Vectors used to store the energies of the graph at its various states.
    std::vector<double> * energyStore = new std::vector<double>;
    std::vector<double> * energyStore2 = new std::vector<double>;
    std::vector<double> * energyStore3 = new std::vector<double>;
    //vectors used to store the dimensionality of the graph at its varius states.
    std::vector<double> * dimensions1 = new std::vector<double>;
    std::vector<double> * dimensions2 = new std::vector<double>;
    std::vector<double> * dimensions3 = new std::vector<double>;

    //sets parameters for the simulation.
    std::cout << "How large of a graph would you like to use? ";
    std::cin >> SIZE;
    std::cout << std::endl << "Input the inverse temperature: ";
    std::cin >> TINV;
    std::cout << "How many sweeps would you like to perform? ";
    std::cin >> maxSweeps;
        
    hGraph * graph = new hGraph(SIZE);
    (*graph) = randomGraph(SIZE);
    //Two graphs are created. One is a complete graph, the other is a random graph with fill probability between 0.25 and 0.75;
    hGraph * graphK = new hGraph(SIZE);
    *graphK = compGraph(SIZE);
    
    hGraph * graphEmpty = new hGraph(SIZE);
    (*graphEmpty) = zeroGraph(SIZE);
    
    std::cout << "Displaying progress for complete graph seed: " << std::endl;
    std::thread threadA(monteCarlo, graph, energyStore, dimensions1, false, "Random Graph");    //Creates a thread to simulate on the random graph
    std::thread threadB(monteCarlo, graphEmpty, energyStore2, dimensions2, false, "Empty Graph");
    monteCarlo(graphK, energyStore3, dimensions3, true, "Complete Graph");                       //Main thread handles the complete graph

    threadA.join();
    threadB.join();
    
    std::cout << "Calculating correlation functions..."  << std::endl;
    
    std::vector<double> * corr1 = new std::vector<double>;
    std::vector<double> * corr2 = new std::vector<double>;
    std::vector<double> * corr3 = new std::vector<double>;
    
    std::thread threadC(correlationFn, energyStore, corr1);
    std::thread threadD(correlationFn, energyStore2, corr2);
    correlationFn(energyStore3, corr3);
    threadC.join();
    threadD.join();

    
    int points = energyStore->size();
    int points2 = energyStore2->size();
    int points3 = energyStore3->size();
    
    
    std::cout << "Graphing energy:" << std::endl;
    
    std::vector<std::vector<double>> xVals;
    std::vector<std::vector<double>> yVals;
    std::vector<double> temp;
    for(int i = 0; i < points; i++)
        temp.push_back(i);
    xVals.push_back(temp);
    yVals.push_back(*energyStore);
    temp.clear();
    for(int i = 0; i < points2; i++)
        temp.push_back(i);
    xVals.push_back(temp);
    yVals.push_back(*energyStore2);
    temp.clear();
    for(int i = 0; i < points3; i++)
        temp.push_back(i);
    
    xVals.push_back(temp);
    yVals.push_back(*energyStore3);
    
    std::cin.clear();
    std::cin.ignore(100, '\n');

    
    //drawMultiGraph(xVals, yVals);
    

    
    yVals.clear();
    yVals.push_back(*corr1);
    yVals.push_back(*corr2);
    yVals.push_back(*corr3);
    
    std::cout << "Graphing correlation function: " << std::endl;
    for (int i = 0; i < 3; i++) {
        for(int j = 0; j < points; j++) {
            corrOut << yVals[i][j];
            corrOut << ',';
            
        }
        corrOut << '\n';
        
    }
    
    //drawMultiGraph(xVals, yVals);

    
    yVals.clear();
    yVals.push_back(*dimensions1);
    yVals.push_back(*dimensions2);
    yVals.push_back(*dimensions3);
    
    std::cout << "Graphing dimensionality: " << std::endl;
   // drawMultiGraph(xVals, yVals);
    
    //deletes pointers
    delete energyStore;
    delete energyStore2;
    delete energyStore3;
    delete corr1;
    delete corr2;
    delete corr3;
}

void monteCarlo (hGraph * graph, std::vector<double> * energy, std::vector<double> * dimensions, bool progress, std::string descriptor) {
    graph->setThreads(NUM_CORES/3);
    std::random_device rd; //c++11 random number generator
    std::mt19937 gen(rd());
    std::uniform_int_distribution<int> randNode(0, SIZE - 1);
    std::uniform_real_distribution<> probDis(0.0, 1.0);
    std::uniform_int_distribution<int> randNum(1, ceil(double(SIZE)/2));
    
    //Parameters needed for individual simulation steps.
    bool run = true;
    int nodeA = 0;
    int nodeB = 0;
    int sweeps = 0;
    int sweepNUM = ((SIZE)*(SIZE-1))/2;//Number of edge flips in a single sweep;
    int n = 0;
    int increases = 0;  //number of times an increase in energy is accepted
    int swaps = 0;      //number of times an edge flip is accepted.
    
    
    simFunction(*graph); //calculates inital graph energy;
    energy->push_back(graph->getHam());

    
    while(run) {
        std::vector<int> xVals;
        std::vector<int> yVals;
        
        int numFlips = randNum(gen);
        int numAccept = 0;
        
        bool isFull = false;
        while(!isFull) {
            while(nodeA == nodeB) { //Randomly selects two nodes;
                nodeA = randNode(gen);
                nodeB = randNode(gen);
            }
            bool same = false;
            for(int l = 0; l < xVals.size(); l++) {
                if(nodeA == xVals[l]) {
                    if(nodeB == yVals[l]) {
                        same = true;
                    }
                }
                
                else if (nodeA == yVals[l]) {
                    if(nodeB == xVals[l]) {
                        same = true;
                    }
                }
            }
            
            if(!same) {
                xVals.push_back(nodeA);
                yVals.push_back(nodeB);
                numAccept++;
            }
            if(numAccept == numFlips) {
                isFull = true;
            }
            nodeA = 0;
            nodeB = 0;
        }
        

        int multiSwap = 0;
        double hamInitial = graph->getHam();
        double hamDiff = simPartial(*graph, xVals, yVals); //gets the change in energy if edge between nodeA and nodeB is flipped
        if(hamDiff <= 0) {  //Accepts change if it decreases energy.
            graph->flipEdge(xVals,yVals);
            graph->acceptPartial(hamDiff);
            swaps++;
            if(xVals.size() > 1) {
                multiSwap++;
            }
        }
        else { //If swap increases energy, the change is accepted given a particular probability.
            
            double prob = exp(TINV*(-hamDiff));
            if(probDis(gen) <= prob) {
                graph->flipEdge(xVals,yVals);
                graph->acceptPartial(hamDiff);
                increases++;
                swaps++;
                if(xVals.size() > 1) {
                    multiSwap++;
                }

            }
            
            
        }
        xVals.clear();
        yVals.clear();
        nodeA = 0;
        nodeB = 0;  //resets
        n += numFlips;        //counts number of attempts
        if(n >= sweepNUM) { //Counts number of sweeps that have been made.
            n = 0;
            sweeps++;
            energy->push_back(graph->getHam());
            graph->setThreads(NUM_CORES/3);
            graph->calcDimension();
            dimensions->push_back(graph->getDimension());
            if((progress == true) && (((sweeps % (maxSweeps/100)) == 0))) {
                std::cout << sweeps << " sweeps completed." << std::endl;
            }

            
        }
        
        
        if(sweeps == maxSweeps) {   //terminates simulation when the maximum number of sweeps has been reached.
            run = false;
        }
        
        
        
    }
    descriptor += ".csv";
    std::ofstream output(descriptor);
    output << *graph;
    output.close();
    std::cout << std::endl;
    std::cout << "Simulation on graph " << descriptor << " complete" << std::endl;

    
    
    
}

void correlationFn(std::vector<double> * data, std::vector <double> * output) {
                                        //Calculates the autocorrelation function when passed a vector with data and a vector to place results
                                        //(pointers). See Monte Carlo Methods in Statistical Physics - 3.21
    int tMax = data->size();
    for(int t = 0; t < tMax; t++ ) {    //The algorithm uses the same variable names (with tp for t') as the algorithm in the book does.
        double sum = 0;                 //Calculates the function at all times. 
        double partialSum = 0;
        double partialSum2 = 0;
        for(int tp = 0; tp < tMax - t; tp++) {
            partialSum += (*data)[tp]*(*data)[tp + t];
        }
        partialSum /= (tMax -t);
        sum += partialSum;
        partialSum = 0;
        
        for(int tp = 0; tp < tMax - t; tp++) {
            partialSum += (*data)[tp];
        }
        partialSum /= (tMax - t);
        
        for(int tp = 0; tp < tMax - t; tp++) {
            partialSum2 += (*data)[tp +t];
        }
        partialSum2 /= (tMax -t);
        
        partialSum *= partialSum2;
        sum -= partialSum;
        if(t > 0) {
            sum /= (*output)[0];
        }
        
        output->push_back(sum);
        
    }
    (*output)[0] = 1;
    
}

