#include <iostream>
#include <string>
#include "hamiltonians.h"
#include "graphics/graphingUtil.hpp"
#include <cmath>
#include <vector>
#include <thread>

void monteCarlo(hGraph * graph, std::vector<double> * energy, std::vector<double> * dimensions, bool progress, std::string descriptor);

int NUM_CORES = 4; //Number of cores in the CPU
double TINV;    //Beta
int SIZE;       //Number of nodes in the graphs
int maxSweeps;  //Number of sweeps that will be performed
int sweepCollect;


//------Simulation Functions------//
std::function<void(hGraph&)> simFunction;
std::function<double(hGraph&,std::vector<int>,std::vector<int>)> simPartial;
//All these variables are declared to be global so that they can be used by multiple threads.



int main() {
    std::ofstream corrOut("correlation.csv");
    std::ofstream dimenOut("dimensionality.csv");
    std::ofstream engOut("energy.csv");
    std::ofstream paramOut("parameters.txt");
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
    std::cout << "Number of available threads? ";
    std::cin >> NUM_CORES;
    
    paramOut << "Graph size: " << SIZE << std::endl;
    paramOut << "Inverse Temperature: " << TINV << std::endl;
    paramOut << "Sweeps performed: " << maxSweeps << std::endl;
    paramOut << "Source term: 0.1" << std::endl;
    paramOut.close();
        
    hGraph * random1 = new hGraph(SIZE); //Three graphs are used, two random ones with a fill fraction between .25 and .75, the other is the empty graph
    (*random1) = randomGraph(SIZE);
    hGraph * random2 = new hGraph(SIZE);
    *random2 = randomGraph(SIZE);
    
    hGraph * graphEmpty = new hGraph(SIZE);
    (*graphEmpty) = zeroGraph(SIZE);
    
    monteCarlo(random1, energyStore, dimensions1, true, "Random Graph 1");    //As the only truly multithreaded calculation is the dimensionality, it is a more efficient use of
    monteCarlo(random2, energyStore2, dimensions2, true, "Random Graph 2");   //processing power to do the individual simulations in series to avoid threads being underutilized.
    monteCarlo(graphEmpty, energyStore3, dimensions3, true, "Empty Graph");

    
    std::cout << "Calculating correlation functions..."  << std::endl;
    
    std::vector<double> * corr1 = new std::vector<double>;
    std::vector<double> * corr2 = new std::vector<double>;
    std::vector<double> * corr3 = new std::vector<double>;
    
    std::thread threadC(correlationFn, dimensions1, corr1);
    std::thread threadD(correlationFn, dimensions2, corr2);
    correlationFn(dimensions3, corr3);
    threadC.join();
    threadD.join();

    
    int points = energyStore->size();
    
    
    std::vector<std::vector<double>> xVals; //These next several lines of code save all the data to CSV files in case there is an issue while plotting.
    std::vector<std::vector<double>> yVals;
    std::vector<double> temp;
    for(int i = 0; i < points; i++)
        temp.push_back(i);
    xVals.push_back(temp);
    yVals.push_back(*energyStore);
    xVals.push_back(temp);
    yVals.push_back(*energyStore2);
    xVals.push_back(temp);
    yVals.push_back(*energyStore3);
    
    for (int i = 0; i < 3; i++) {
        for(int j = 0; j < points; j++) {
            engOut << yVals[i][j];
            if(j != points-1) {
                engOut << ',';
            }
            
        }
        engOut << '\n';
        
    }
    engOut.close();
    
    yVals.clear();
    yVals.push_back(*corr1);
    yVals.push_back(*corr2);
    yVals.push_back(*corr3);
    
    for (int i = 0; i < 3; i++) {
        for(int j = 0; j < points; j++) {
            corrOut << yVals[i][j];
            if(j != points-1) {
                corrOut << ',';
            }
            
        }
        corrOut << '\n';
        
    }
    corrOut.close();

    
    yVals.clear();
    yVals.push_back(*dimensions1);
    yVals.push_back(*dimensions2);
    yVals.push_back(*dimensions3);
    
    for (int i = 0; i < 3; i++) {
        for(int j = 0; j < points; j++) {
            dimenOut << yVals[i][j];
            if(j != points-1) {
                dimenOut << ',';
            }
        }
        dimenOut << '\n';
        
    }
    dimenOut.close();

    
    std::cin.clear();
    std::cin.ignore(100, '\n');
    
    yVals.clear();
    
    yVals.push_back(*energyStore);
    yVals.push_back(*energyStore2);
    yVals.push_back(*energyStore3);
    
    std::cout << "Graphing energy..." << std::endl; //plots energy
    
    drawMultiGraph(xVals, yVals);
    

    
    yVals.clear();
    yVals.push_back(*corr1);
    yVals.push_back(*corr2);
    yVals.push_back(*corr3);
    
    std::cout << "Graphing correlation function... " << std::endl; //plots correlation function
    drawMultiGraph(xVals, yVals);

    
    yVals.clear();
    yVals.push_back(*dimensions1);
    yVals.push_back(*dimensions2);
    yVals.push_back(*dimensions3);
    
    
    std::cout << "Graphing dimensionality... " << std::endl; //plots dimensionality.
    drawMultiGraph(xVals, yVals);
    
    std::cout << "Be sure to move all result files to a different folder to ensure they are not overwritten." << std::endl;
    
    //deletes pointers
    delete energyStore;
    delete energyStore2;
    delete energyStore3;
    delete corr1;
    delete corr2;
    delete corr3;
}

void monteCarlo (hGraph * graph, std::vector<double> * energy, std::vector<double> * dimensions, bool progress, std::string descriptor) {
    graph->setThreads(NUM_CORES);
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
    
    graph->calcDimension();
    dimensions->push_back(graph->getDimension());

    
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
            graph->setThreads(NUM_CORES);
            graph->calcDimension();
            dimensions->push_back(graph->getDimension());
            if((maxSweeps >= 100 ) && (progress == true) && (((sweeps % (maxSweeps/100)) == 0))) {
                std::cout << descriptor << ": " << sweeps << " sweeps completed." << std::endl;
            }

            
        }
        
        
        if(sweeps == maxSweeps) {   //terminates simulation when the maximum number of sweeps has been reached.
            run = false;
        }
        
        
        
    }
    std::cout << "Simulation on graph " << descriptor << " complete" << std::endl;
    descriptor += ".csv";
    std::ofstream output(descriptor);
    output << *graph;
    output.close();
    std::cout << std::endl;

    
    
    
}


