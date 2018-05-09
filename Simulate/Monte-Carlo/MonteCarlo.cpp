#include <iostream>
#include <string>
#include "hamiltonians.h"
#include <cmath>
#include <vector>
#include "TGraphErrors.h"
#include "TCanvas.h"
#include "TMultiGraph.h"
#include <thread>

void monteCarlo(hGraph * graph, std::vector<double> * energy);
void correlationFn(std::vector<double> * data, std::vector <double> * output); //Simulation itself must be a function in order to be paralelized.


double TINV;    //Beta
int SIZE;       //Number of nodes in the graphs
int maxSweeps;  //Number of sweeps that will be performed


//------Simulation Functions------//
std::function<void(hGraph&)> simFunction;
std::function<double(hGraph&,int, int)> simPartial;
//All these variables are declared to be global so that they can be used by multiple threads.



int main() {
    //Sets simulation functions
    simFunction = basicSquareHam;
    simPartial = basicSquarePartial;

    //Vectors used to store the energies of the graph at its various states.
    std::vector<double> * energyStore = new std::vector<double>;
    std::vector<double> * energyStore2 = new std::vector<double>;
    std::vector<double> * corrFn = new std::vector<double>;
    
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
    
    
    std::thread threadA(monteCarlo, graph, energyStore);    //Creates a thread to simulate on the random graph
    monteCarlo(graphK, energyStore2);                       //Main thread handles the complete graph
    
    threadA.join();
    
    int points = energyStore->size();
    int points2 = energyStore2->size();

    Double_t x[points], y[points], x2[points2], y2[points2], corrY[points]; //Declares arrays for use as data values when drawing plots.
    correlationFn(energyStore, corrFn); //calculates the correlation function for all time steps


    for(int i = 0; i < points; i++) {
        x[i] = i;
        y[i] = (*energyStore)[i];   //fills the data arrays for use when drawing plots. Note, this is necessary because the energy storage containers are vectors,
                                    //but Root takes arrays as parameters when creating gaphs
    }
    
    for(int i = 0; i < points2; i++) {
        x2[i] = i;
        y2[i] = (*energyStore2)[i];
    }
    for(int i = 0; i < points; i++) {
        corrY[i] = (*corrFn)[i];
        
    }
    
    
    //-----------------------------Graph Drawing-----------------------//
    TCanvas *c1 = new TCanvas("c1","Graph Examples",200,10,600,400);

    //Specifics here are not important
    TMultiGraph *graphMulti = new TMultiGraph();
    graphMulti->SetTitle("Energy"); //Set the title to be drawn on the graph
    
    
    TGraph *gr1 = new TGraph (points, x, y);
    gr1->SetLineColor(2);
    gr1->SetLineWidth(1);
    gr1->SetFillStyle(3005);
    //Sets display preferences for the two data series.
    TGraph *gr2 = new TGraph (points, x2, y2);
    gr2->SetLineColor(4);
    gr2->SetLineWidth(1);
    gr2->SetFillStyle(3005);

    
    graphMulti->Add(gr1);
    graphMulti->Add(gr2);
    //Adds the data series to the window.
    graphMulti->Draw("AC");
    c1->Print("test.png");
    //draws the graph and outputs it to a file.
    
    //The following does everything above but for the correlation function.
    TCanvas *c2 = new TCanvas("c2","Autocorrelation Examples",200,10,600,400);
    
    TGraph * corrGraph = new TGraph(points, x, corrY);
    corrGraph->SetTitle("Autocorrelation function");
    
    corrGraph->SetLineColor(2);
    corrGraph->SetLineWidth(1);
    corrGraph->SetFillStyle(3005);

    corrGraph->Draw("AC");
    c2->Print("correlation.png");

    //deletes pointers
    delete energyStore;
    delete energyStore2;
    delete corrFn;
    
}

void monteCarlo (hGraph * graph, std::vector<double> * energy ) {
    std::random_device rd; //c++11 random number generator
    std::mt19937 gen(rd());
    std::uniform_int_distribution<int> randNode(0, SIZE - 1);
    std::uniform_real_distribution<> probDis(0.0, 1.0);
    //Parameters needed for individual simulation steps.
    bool run = true;
    int nodeA = 0;
    int nodeB = 0;
    int sweeps = 0;
    int sweepNUM = ((SIZE)*(SIZE-1))/2;//Number of edge flips in a single sweep;
    int n = 0;
    int increases = 0;  //number of times an increase in energy is accepted
    int swaps = 0;      //number of times an edge flip is accepted.
    

    
    while(run) {
        
        while(nodeA == nodeB) { //Randomly selects two nodes;
            nodeA = randNode(gen);
            nodeB = randNode(gen);
            
        }
        
        simFunction(*graph); //calculates inital graph energy;

        double hamInitial = graph->getHam();
        double hamDiff = simPartial(*graph, nodeA, nodeB); //gets the change in energy if edge between nodeA and nodeB is flipped
        if(hamDiff <= 0) {  //Accepts change if it decreases energy.
            graph->flipEdge(nodeA,nodeB);
            graph->acceptPartial(hamDiff);
            swaps++;
        }
        else { //If swap increases energy, the change is accepted given a particular probability.
            
            double prob = exp(TINV*(-hamDiff));
            if(probDis(gen) <= prob) {
                graph->flipEdge(nodeA,nodeB);
                graph->acceptPartial(hamDiff);
                increases++;
                swaps++;

            }
            
            
        }
        
        nodeA = 0;
        nodeB = 0;  //resets
        n++;        //counts number of attempts
        if(n == sweepNUM) { //Counts number of sweeps that have been made.
            n = 0;
            sweeps++;
            energy->push_back(graph->getHam());
            if((maxSweeps >= 10) && sweeps % (maxSweeps/10) == 0) {
                std::cout << sweeps << " sweeps performed" << std::endl;
            }
            
        }
        
        if(sweeps == maxSweeps) {   //terminates simulation when the maximum number of sweeps has been reached.
            run = false;
        }
        
        
        
    }
    std::cout << std::endl;
    std::cout << *graph << std::endl;
    std::cout << "Moves that increased energy: " << increases << std::endl;
    std::cout << "total swaps accepted: " << swaps <<std::endl;

    
    
    
}

void correlationFn(std::vector<double> * data, std::vector <double> * output) {
    std::cout << "Calculating autocorrelation function..." << std::endl;    //Calculates the autocorrelation function when passed a vector with data and a vector to place results
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
        
        output->push_back(sum);
        
    }
    
}

