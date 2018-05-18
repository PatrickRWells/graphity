#include <iostream>
#include <string>
#include "hamiltonians.h"
#include "graphics/graphUtil/graphingUtil.hpp"
#include "graphics/graphImager/graphImager.h"
#include <cmath>
#include <vector>
#include <thread>

void monteCarlo(hGraph * graph, std::vector<bool> observe, std::vector<double> data[][NUM_OBSERVABLES], int simNum, bool progress, std::string descriptor);

int NUM_CORES = 4; //Number of cores in the CPU
double TINV;    //Beta
int SIZE;       //Number of nodes in the graphs
int maxSweeps;  //Number of sweeps that will be performed
int sweepCollect;

const int DIMEN = 0;
const int DIMEN_CORR = 1; //These constants determine which vector from the array declared below to look at when looking for a particular data set
const int ENERGY_CORR = 2;
const int USER_IN = 3;
const int ENERGY = 4;
const int AVG_DEGREE = 5;


//------Simulation Functions------//
std::function<void(hGraph&)> simFunction;
std::function<double(hGraph&,std::vector<int>,std::vector<int>)> simPartial;
//All these variables are declared to be global so that they can be used by multiple threads.



int main() {
    
    std::vector<bool> observables (NUM_OBSERVABLES, false); //true/false. Tells the program which data to save.
    
    std::ofstream engOut("energy.csv");
    std::ofstream paramOut("parameters.txt");

    std::ofstream dimenCorrOut;
    std::ofstream dimenOut; //Output files streams (may or may not be used depending on what the user actuall needs);
    std::ofstream engCorrOut;
    std::ofstream aveDegOut;
    std::ofstream allOut;
    
    
    //Sets simulation functions
    simFunction = basicSquareHam;
    simPartial = basicSquarePartial;


    while(true) {
        std::cout << "Would you like to input a graph from a file? (y/n) "; //Checks if the user would like to input a graph from a file (vs. using randomly initialized graphs);
        char c = getchar();
        if (toupper(c) == 'Y') {
            observables[USER_IN] = true;
            std::cin.clear();
            std::cin.ignore(100, '\n');
            
            break;
        }
        else if(toupper(c) == 'N') {
            std::cin.clear();
            std::cin.ignore(100, '\n');
            
            break;
        }
        else {
            std::cout << "Invalid input." << std::endl;
            std::cin.clear();
            std::cin.ignore(100, '\n');
            
        }
    }
    
    hGraph ** graphs;
    int numGraphs = 0;
    
    if(observables[USER_IN]) { //reads in the graph file if the user chooses to use one.
        hGraph * temp = readGraphFile(numGraphs);
        graphs = new hGraph *[numGraphs];
        for(int i = 0; i < numGraphs; i++) {
            graphs[i] = new hGraph;
            *graphs[i] =temp[i];
        }
       
        SIZE = graphs[0]->getSize();
        std::cout << "Graph size detected to be " << SIZE << std::endl;
    }
    else {
        numGraphs = 3;
    }
    std::vector<double> data[numGraphs][NUM_OBSERVABLES];
    
    //sets parameters for the simulation.
    if(!observables[USER_IN]) {
    std::cout << "How large of a graph would you like to use? ";
    std::cin >> SIZE;
    }
    std::cout << "Input the inverse temperature: ";
    std::cin >> TINV;
    std::cout << "How many sweeps would you like to perform? ";
    std::cin >> maxSweeps;
    std::cout << "Number of available threads? ";
    std::cin >> NUM_CORES;
    
    std::cin.clear();
    std::cin.ignore(100, '\n');

    
    while(true) { //Prompts the user regarding which quantities should be calculated (and therefore plotted as well);
        std::cout << "Would you like to calculate the dimensionality? (y/n) ";
        char c = getchar();
        if (toupper(c) == 'Y') {
            observables[DIMEN] = true;
            std::cin.clear();
            std::cin.ignore(100, '\n');

            break;
        }
        else if(toupper(c) == 'N') {
            std::cin.clear();
            std::cin.ignore(100, '\n');

            break;
        }
        else {
            std::cout << "Invalid input." << std::endl;
            std::cin.clear();
            std::cin.ignore(100, '\n');

        }
    }
    if(observables[DIMEN]) {
        while(true) {
            std::cout << "Would you like to calculate the correlation function for the dimensionality? (y/n) ";
            char c = getchar();
            if (toupper(c) == 'Y') {
                observables[DIMEN_CORR] = true;
                std::cin.clear();
                std::cin.ignore(100, '\n');
                
                break;
            }
            else if(toupper(c) == 'N') {
                std::cin.clear();
                std::cin.ignore(100, '\n');
                
                break;
            }
            else {
                std::cout << "Invalid input." << std::endl;
                std::cin.clear();
                std::cin.ignore(100, '\n');
                
            }
        }
    }
    while(true) {
        std::cout << "Would you like to calculate the correlation function for the energy? (y/n) ";
        char c = getchar();
        if (toupper(c) == 'Y') {
            observables[ENERGY_CORR] = true;
            std::cin.clear();
            std::cin.ignore(100, '\n');
            
            break;
        }
        else if(toupper(c) == 'N') {
            std::cin.clear();
            std::cin.ignore(100, '\n');
            
            break;
        }
        else {
            std::cout << "Invalid input." << std::endl;
            std::cin.clear();
            std::cin.ignore(100, '\n');
            
        }
    }
    
    while(true) {
        std::cout << "Would you like to calculate the average node degree? (y/n) ";
        char c = getchar();
        if (toupper(c) == 'Y') {
            observables[AVG_DEGREE] = true;
            std::cin.clear();
            std::cin.ignore(100, '\n');
            
            break;
        }
        else if(toupper(c) == 'N') {
            std::cin.clear();
            std::cin.ignore(100, '\n');
            
            break;
        }
        else {
            std::cout << "Invalid input." << std::endl;
            std::cin.clear();
            std::cin.ignore(100, '\n');
            
        }
    }
    
    
    
    paramOut << "Graph size: " << SIZE << std::endl;
    paramOut << "Inverse Temperature: " << TINV << std::endl;
    paramOut << "Sweeps performed: " << maxSweeps << std::endl;
    paramOut << "Source term: -0.1" << std::endl;
    paramOut.close(); //This parameter file contians all the simulation's parameters for future reference.
    
    

    std::cin.clear();
    
    
    std::string description[numGraphs];
    if(observables[USER_IN]) {
        for(int i = 0; i < numGraphs; i++){
            std::cout << "Enter a short descriptor for graph " << i << ": "; //Prompts the users for ways to describe the graph if the graphs were input from a file.
            getline(std::cin, description[i]);
        }

        for(int i = 0; i < numGraphs; i++) { //Runs the monte-carlo simulations
            monteCarlo(graphs[i], observables, data, i, true, description[i]);
            
        }
        
    }

    else { //If the user did not input the graphs from a file, two random graphs and an empty graph are used as seeds.
        graphs = new hGraph * [3];
        graphs[0] = new hGraph(SIZE);
        *graphs[0] = randomGraph(SIZE);
        graphs[1] = new hGraph(SIZE);
        *graphs[1] = randomGraph(SIZE);
        graphs[2] = new hGraph(SIZE);
        *graphs[2] = zeroGraph(SIZE);
        numGraphs = 3;
        description[0] = "Random Graph 1";
        description[1] = "Random Graph 2";
        description[2] = "Empty Graph";

        monteCarlo(graphs[0], observables, data, 0, true, description[0]);    //As the only truly multithreaded calculation is the dimensionality, it is a more efficient use of
        monteCarlo(graphs[1], observables, data, 1, true, description[1]);   //processing power to do the individual simulations in series to avoid threads being underutilized.
        monteCarlo(graphs[2], observables, data, 2, true, description[2]);
    }
    
    
    
    
    
    allOut.open("allOut.csv");
    allOut << SIZE << std::endl;
    for(int i = 0; i < numGraphs; i++) {
        allOut << *(graphs[0]); //Outputs all the graphs to a single file. Graphs are also outputted to individual files in the monte-carlo simulation.
    }
    allOut.close();
    
    
    if(observables[DIMEN_CORR]) {
        std::cout << "Calculating dimensionality corrrelation function..."  << std::endl;
        for(int i = 0; i < numGraphs; i++) {
            correlationFn(data, i, DIMEN, DIMEN_CORR);
        }
        
        std::cout << "Dimensionality correlation function calculation complete" << std::endl;
    }
    
    
    
    if(observables[ENERGY_CORR]) {
        std::cout << "Calculating energy correlation function..." << std::endl;
            for(int i = 0; i < numGraphs; i++) {
                correlationFn(data, i, ENERGY, ENERGY_CORR);
            }
        std::cout << "Energy correlation function calculation complete" << std::endl;
    }
        
        
        
    
    
    for (int i = 0; i < numGraphs; i++) { //Outputs energy to a CSV file
        for(int j = 0; j < data[i][ENERGY].size(); j++) {
            engOut << data[i][ENERGY][j];
            if(j != data[i][ENERGY].size()-1) {
                engOut << ',';
            }
            
        }
        engOut << '\n';
        
    }
    engOut.close();
    
    if(observables[DIMEN_CORR]) { //Outputs dimensionality correlation function to a CSV file
        dimenCorrOut.open("dimenCorrelation.csv");
        for (int i = 0; i < numGraphs; i++) {
            for(int j = 0; j < data[i][DIMEN_CORR].size(); j++) {
                dimenCorrOut << data[i][DIMEN_CORR][j];
                if(j != data[i][DIMEN_CORR].size() - 1) {
                    dimenCorrOut << ',';
                }
                
            }
            dimenCorrOut << '\n';
        }
        dimenCorrOut.close();
    }

    if(observables[DIMEN]) { //outputs dimensionality to a CSV file
        dimenOut.open("dimensionality.csv");
        for (int i = 0; i < numGraphs; i++) {
            for(int j = 0; j < data[i][DIMEN].size(); j++) {
                dimenOut << data[i][DIMEN][j];
                if(j != data[i][DIMEN].size()-1) {
                    dimenOut << ',';
                }
            }
            dimenOut << '\n';
        }
        dimenOut.close();
    }
    
    if(observables[ENERGY_CORR]) { //Outputs the energy correlation function to a CSV file
        engCorrOut.open("energyCorrelation.csv");
        for (int i = 0; i < numGraphs; i++) {
            for(int j = 0; j < data[i][ENERGY_CORR].size(); j++) {
                engCorrOut << data[i][ENERGY_CORR][j];
                if(j != data[i][ENERGY_CORR].size()-1) {
                    engCorrOut << ',';
                }
            }
            engCorrOut << '\n';
        }
        engCorrOut.close();
    }
    
    
    
    
    std::cout << "Data outputed to CSV files" << std::endl;

    
    std::cin.clear();
    //std::cin.ignore(100, '\n');
    
    
    std::cout << "Graphing energy..." << std::endl; //plots energy
    drawMultiGraph(data, numGraphs, ENERGY);
    
    
    if(observables[ENERGY_CORR]) {
        std::cout << "Graphing energy correlation function... " << std::endl; //plots energy correlation function
        drawMultiGraph(data, numGraphs, ENERGY_CORR);
        
    }

    if(observables[DIMEN]) {
        std::cout << "Graphing dimensionality... " << std::endl; //plots dimensionality.
        drawMultiGraph(data, numGraphs, DIMEN);
    }
    
    if(observables[DIMEN_CORR]) {
        std::cout << "Graph dimensionality correlation function... " << std::endl; //plots dimensionality correlation function
        drawMultiGraph(data, numGraphs, DIMEN_CORR);
    }
    
    if(observables[AVG_DEGREE]) {
        std::cout << "Graph average node degree" << std::endl; //plots dimensionality correlation function
        drawMultiGraph(data, numGraphs, AVG_DEGREE);
        
    }
    
    while(true) {
        std::cout << "Would you like to make images of the graphs? (y/n) "; //Checks if the user would like to input a graph from a file (vs. using randomly initialized graphs);
        char c = getchar();
        if (toupper(c) == 'Y') {
            std::cin.clear();
            std::cin.ignore(100, '\n');
            for(int i = 0; i < numGraphs; i++) {
                std::cout << "Imaging graph " << description[i] << "..." << std::endl;
                graphImage(*graphs[i]);
                
            }
            
            break;
        }
        else if(toupper(c) == 'N') {
            std::cin.clear();
            std::cin.ignore(100, '\n');
            
            break;
        }
        else {
            std::cout << "Invalid input." << std::endl;
            std::cin.clear();
            std::cin.ignore(100, '\n');
            
        }
    }

    
    
    std::cout << "Be sure to move all result files to a different folder to ensure they are not overwritten." << std::endl;
    
}

void monteCarlo (hGraph * graph, std::vector<bool> observe, std::vector<double> data[][NUM_OBSERVABLES], int simNum, bool progress, std::string descriptor) {
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
    data[simNum][ENERGY].push_back(graph->getHam());

    if(observe[DIMEN]) {
        graph->calcDimension();
        data[simNum][DIMEN].push_back(graph->getDimension());
    }
    
    if(observe[AVG_DEGREE]) {
        double sum = 0;
        for(int i = 0; i < SIZE; i++) {
            sum += graph->getDegree(i);
        }
        data[simNum][AVG_DEGREE].push_back(sum/SIZE);
    }

    
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
            data[simNum][ENERGY].push_back(graph->getHam());
            
            if(observe[DIMEN]) {
                graph->setThreads(NUM_CORES);
                graph->calcDimension();
                data[simNum][DIMEN].push_back(graph->getDimension());
                
            }
            
            if(observe[AVG_DEGREE]) {
                double sum = 0;
                for(int i = 0; i < SIZE; i++) {
                    sum += graph->getDegree(i);
                }
                data[simNum][AVG_DEGREE].push_back(sum/SIZE);
            }
            
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
    output << graph->getSize() << std::endl;
    output << *graph;
    output.close();
    std::cout << std::endl;

    
    
    
}


