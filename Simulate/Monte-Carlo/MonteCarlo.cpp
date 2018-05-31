#include <iostream>
#include <string>
#include "hamiltonians.h"
#include "graphics/graphUtil/graphingUtil.hpp"
#include "graphics/graphImager/graphImager.h"
#include <cmath>
#include <vector>
#include <thread>
#include <mutex>

void monteCarlo(hGraph * graph, std::vector<bool> observe, std::vector<double> ** data, int simNum, bool progress, std::string descriptor);
void wolffAlgorithm(hGraph * graph, std::vector<bool> observe, std::vector<double> data[][NUM_OBSERVABLES], int simNum, bool progress, std::string descriptor);
std::mutex mut;


bool getTF();

bool WOLFF = false;
int NUM_CORES = 4; //Number of cores in the CPU
double TINV;    //Beta
int SIZE;       //Number of nodes in the graphs
int maxSweeps;  //Number of sweeps that will be performed
int collectionTime = 1;
bool correlationCollect = false;




//------Simulation Functions------//
std::function<void(hGraph&)> simFunction;
std::function<double(hGraph&,std::vector<int>,std::vector<int>)> simPartial;
//All these variables are declared to be global so that they can be used by multiple threads.


int main() {
    
    
    
    std::vector<bool> observables (NUM_OBSERVABLES, false); // true/false. Tells the program which data to save.
    std::vector<bool> plot(NUM_OBSERVABLES, false);
    
    std::ofstream engOut("energy.csv");
    std::ofstream paramOut("parameters.txt");

    std::ofstream dimenCorrOut;
    std::ofstream dimenOut; //Output files streams (may or may not be used depending on what the user actuall needs);
    std::ofstream engCorrOut;
    std::ofstream avgDegOut;
    std::ofstream allOut;
    std::ofstream avgDegCorrOut;
    std::ofstream eulerCharOut;
    
    //Sets simulation functions
    simFunction = basicSquareHam;
    simPartial = basicSquarePartial;


    std::cout << "Would you like to input a graph from a file? (y/n) "; //Checks if the user would like to input a graph from a file (vs. using randomly initialized graphs);
    observables[USER_IN] = getTF(); //The constants corresponding to parameters (such as USER_IN) are defined in graphingUtil.hpp

    
    hGraph ** graphs;
    int numGraphs = 0;
    
    if(observables[USER_IN]) { //reads in the graph file if the user chooses to use one.
        readGraphFile(&graphs, numGraphs);
        SIZE = graphs[0]->getSize();
        std::cout << "Graph size detected to be " << SIZE << std::endl;
        std::cout << "Would you like to manually set how often to collect data? (y/n)";
        correlationCollect = getTF();
        if(correlationCollect) {
            std::string temp;
            std::cout << "Input how often you would like to collect data (in sweeps): ";
            std::cin >> temp;
            collectionTime = std::stoi(temp);
        }
    }
       
    else {
        numGraphs = 3;
    }
    std::vector<double> ** data = new std::vector<double> * [numGraphs];
    for(int i = 0; i < numGraphs; i++) {
        data[i] = new std::vector<double> [NUM_OBSERVABLES];
        
    }
    
    //sets parameters for the simulation.
    if(!observables[USER_IN]) {
    std::cout << "How large of a graph would you like to use? ";
    std::cin >> SIZE;
    }
    std::cout << "Input the source term: ";
    std::cin >> SOURCE;
    std::cout << "Input the inverse temperature: ";
    std::cin >> TINV;
    std::cout << "How many sweeps would you like to perform? ";
    std::cin >> maxSweeps;
    std::cout << "Number of available threads? ";
    std::cin >> NUM_CORES;
    
    std::cin.clear();
    std::cin.ignore(100, '\n');

    
    std::cout << "Would you like to calculate the dimensionality? (y/n) ";
    observables[DIMEN] = getTF();

    if(observables[DIMEN]) {
        std::cout << "Would you like to calculate the correlation function for the dimensionality? (y/n) ";
        observables[DIMEN_CORR] = getTF();
    }

    std::cout << "Would you like to calculate the correlation function for the energy? (y/n) ";
    observables[ENERGY_CORR] = getTF();

    
    std::cout << "Would you like to calculate the average node degree? (y/n) ";
    observables[AVG_DEGREE] = getTF();

    std::cout << "Would you like to calculate the correlation function for the average node degree? (y/n) ";
    observables[AVG_DEGREE_CORR] = getTF();
    
    std::cout << "Would you like to calculate the Euler Characteristic for the graph? (y/n) ";
    observables[EULER_CHAR] = getTF();
    
    
    paramOut << "Graph size: " << SIZE << std::endl;
    paramOut << "Inverse Temperature: " << TINV << std::endl;
    paramOut << "Sweeps performed: " << maxSweeps << std::endl;
    paramOut << "Source term: " << SOURCE << std::endl;
    paramOut.close(); //This file contians all the simulation's parameters for future reference.
    
    

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
        description[0] = "Random Graph 1";
        description[1] = "Random Graph 2";
        description[2] = "Empty Graph";

        
        std::thread threadA(monteCarlo, graphs[0], observables, data, 0, true, description[0]);   //processing power to do the individual simulations in series to avoid threads being underutilized.
        std::thread threadB(monteCarlo, graphs[1], observables, data, 1, true, description[1]);   //processing power to do the individual simulations in series to avoid threads being underutilized.
        monteCarlo(graphs[2], observables, data, 2, true, description[2]);
        threadA.join();
        threadB.join();
    }
    
    
    
    
    
    allOut.open("allOut.csv");
    allOut << SIZE << std::endl;
    for(int i = 0; i < numGraphs; i++) {
        allOut << *(graphs[0]); //Outputs all the graphs to a single file. Graphs are also outputted to individual files in the monte-carlo simulation.
    }
    allOut.close();
    
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
    
    //plots various data sets
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
    
    if(observables[AVG_DEGREE_CORR]) {
        std::cout << "Calculating average degree correlation function..." << std::endl;
        for(int i = 0; i < numGraphs; i++) {
            correlationFn(data, i , AVG_DEGREE, AVG_DEGREE_CORR);
        }
        std::cout << "Average degree correlation function calculation complete" << std::endl;
    }
    
    
    
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
    
    if(observables[AVG_DEGREE]) {
        avgDegOut.open("averageDegree.csv");
        for (int i = 0; i < numGraphs; i++) {
            for(int j = 0; j < data[i][AVG_DEGREE].size(); j++) {
                avgDegOut << data[i][AVG_DEGREE][j];
                if(j != data[i][AVG_DEGREE].size()-1) {
                    avgDegOut << ',';
                }
            }
            avgDegOut << '\n';
        }
        avgDegOut.close();
    }
    
    if(observables[AVG_DEGREE_CORR]) {
        avgDegCorrOut.open("averageDegreeCorrelation.csv");
        for (int i = 0; i < numGraphs; i++) {
            for(int j = 0; j < data[i][AVG_DEGREE_CORR].size(); j++) {
                avgDegCorrOut << data[i][AVG_DEGREE_CORR][j];
                if(j != data[i][AVG_DEGREE_CORR].size()-1) {
                    avgDegCorrOut << ',';
                }
            }
            avgDegCorrOut << '\n';
        }
        avgDegCorrOut.close();
    }

    if(observables[EULER_CHAR]) {
        eulerCharOut.open("eulerChar.csv");
        for (int i = 0; i < numGraphs; i++) {
            for(int j = 0; j < data[i][EULER_CHAR].size(); j++) {
                eulerCharOut << data[i][EULER_CHAR][j];
                if(j != data[i][EULER_CHAR].size()-1) {
                    eulerCharOut << ',';
                }
            }
            eulerCharOut << '\n';
        }
        eulerCharOut.close();
    }
    
    
    
    
    std::cout << "Data outputed to CSV files" << std::endl;
    
    
    std::cin.clear();
    
    std::cout << "Would you like to plot the energy? (y/n) ";
    
    if(getTF()) {
        std::cout << "Graphing energy..." << std::endl; //plots energy
        drawMultiGraph(data, numGraphs, ENERGY);
    }
    
    if(observables[ENERGY_CORR]) {
        std::cout << "Would you like to plot the energy correlation function? (y/n) ";
        if(getTF()) {
            std::cout << "Graphing energy correlation function... " << std::endl; //plots energy correlation function
            drawMultiGraph(data, numGraphs, ENERGY_CORR);
        }
        
    }

    if(observables[DIMEN]) {
        std::cout << "Would you like to plot the dimensionality? (y/n) ";
        if(getTF() ){
            std::cout << "Graphing dimensionality... " << std::endl; //plots dimensionality.
            drawMultiGraph(data, numGraphs, DIMEN);
        }
    }
    
    if(observables[DIMEN_CORR]) {
        std::cout << "Would you like to plot the dimensionality correlation function? (y/n) ";
        if(getTF()) {
            std::cout << "Graph dimensionality correlation function... " << std::endl; //plots dimensionality correlation function
            drawMultiGraph(data, numGraphs, DIMEN_CORR);
        }
    }
    
    if(observables[AVG_DEGREE]) {
        std::cout << "Would you like to plot the average node degree? (y/n) ";
        if(getTF()) {
            std::cout << "Graph average node degree" << std::endl; //
            drawMultiGraph(data, numGraphs, AVG_DEGREE);
            
        }
        
    }
    
    if(observables[AVG_DEGREE_CORR]) {
        std::cout << "Would you like to plot the average node degree correlation function? (y/n) ";
        if(getTF()) {
            std::cout << "Graph average node degree correlation function" << std::endl; //plots dimensionality correlation function
            drawMultiGraph(data, numGraphs, AVG_DEGREE_CORR);
            
        }
        
    }
    
    if(observables[EULER_CHAR]) {
        std::cout << "Would you like to plot the Euler characteristic? (y/n) ";
        if(getTF()) {
            std::cout << "Graphing Euler Characteristic" << std::endl;
            drawMultiGraph(data, numGraphs, EULER_CHAR);
        }
        
        
    }
    
    std::cout << "Would you like to make images of the graphs? (y/n) "; //Asks the user if they want to create a fancy PNG of the graph. See graphImager
    if(getTF()) {
        for(int i = 0; i < numGraphs; i++) {
            std::cout << "Imaging graph " << description[i] << "..." << std::endl;
            graphImage(*graphs[i]);
        }
    }

    
    delete [] data;
    std::cout << "Be sure to move all result files to a different folder to ensure they are not overwritten." << std::endl;
    
}
void wolffAlgorithm(hGraph * graph, std::vector<bool> observe, std::vector<double> data[][NUM_OBSERVABLES], int simNum, bool progress, std::string descriptor) {
    std::cout << *graph << std::endl;
    int selected = 0;
    graph->setThreads(NUM_CORES);
    std::random_device rd; //c++11 random number generator
    std::mt19937 gen(rd());
    std::uniform_int_distribution<int> randNode(0, SIZE - 1);
    std::uniform_real_distribution<> probDis(0.0, 1.0);
    
    bool run = true;
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
        int nodeA = 0;
        int nodeB = 0;

        int nonAccept = 0;
        std::vector<intPair> values;
        std::vector<intPair> attempted;
        bool on;
        while(nodeA == nodeB) {
            nodeA = randNode(gen);
            nodeB = randNode(gen);
        }
        on = graph->isConnected(nodeA, nodeB);
        intPair tempPair(nodeA, nodeB);
        values.push_back(tempPair);
        attempted.push_back(tempPair);
        std::cout << values.size() << std::endl;

        intPair tempA;
        intPair tempB;
        double prob = 1 - exp(-TINV);
        std::cout << prob << std::endl;
        for(int a = 0; a < values.size(); a++) {
            for(int j = 0; j < SIZE; j++) {
                nodeA = values[a].getPair()[0];
                nodeB = values[a].getPair()[1];
                tempA = intPair(nodeA, j);
                tempB = intPair(nodeB, j);
                if(nodeA != j && (graph->isConnected(nodeA, j) == on)) {
                    bool tested = true;
                    for(int t = 0; t < attempted.size(); t++) {
                        if(attempted[t].equals(tempA)) {
                            break;
                        }
                        if(t == (attempted.size() - 1)) {
                            tested = false;
                        }
                        
                    }
                    if(!tested) {
                        attempted.push_back(tempA);
                        if(probDis(gen) <= prob) {
                            values.push_back(tempA);
                        }
                    }
                }
                if(nodeB != j && (graph->isConnected(nodeB, j) == on)) {
                    bool tested = true;
                    for(int t = 0; t < attempted.size(); t++) {
                        if(attempted[t].equals(tempB)) {
                            break;
                        }
                        if(t == (attempted.size() - 1)) {
                            tested = false;
                        }
                        
                    }
                    if(!tested) {
                        attempted.push_back(tempB);
                        if(probDis(gen) <= prob) {
                            values.push_back(tempB);
                        }
                    }
                }
                
                
            }
        }
        
        for(int i = 0; i < attempted.size(); i++) {
            attempted[i].print();
            std::cout << std::endl;
        }
        
        std::cout << std::endl;
        
        for(int i = 0; i < values.size(); i++) {
            values[i].print();
            std::cout << std::endl;
        }
        run = false;
    }
}

void monteCarlo (hGraph * graph, std::vector<bool> observe, std::vector<double> ** data, int simNum, bool progress, std::string descriptor) {
    int selected = 0;
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
        bool select = true;
        hGraph temp(graph->getSize());
        temp = *graph;
        temp.flipEdge(xVals, yVals);
        
        //If two graphs are isomorphic, it's like it never happened.
        //If they are not isomorphic, Calculate the automorphism group sizes. If the size of new is greater than size of old, just select the new graph
        //If the size of old is greater than size of the new, generate a random number. If that number is less than (new group size)/(old group size)

        mut.lock();
        bool iso = isIsomorphic(*graph, temp);
        mut.unlock();

        if (iso) {
            continue;
        }
        
    
        else {
            mut.lock();
            std::vector<double> aGrp1 = graph->autoGroupSize();
            std::vector<double> aGrp2 = temp.autoGroupSize();
            mut.unlock();
            double randAuto = (aGrp2[0]/aGrp1[0])*pow(10, aGrp2[1] - aGrp1[1]);
            //std::cout << aGrp2[0]*pow(10, aGrp2[1]) << ", " << aGrp1[0]*pow(10, aGrp1[1]) << ", " << randAuto << std::endl;

            
            if(randAuto >= 1.0) {
                selected++;
                select = true;
            }
            else {
                
                //std::cout << randAuto << std::endl;
                if(probDis(gen) <= randAuto) {
                    selected++;
                    select = true;
                }
            }
            
        }
        

        if(select) {
            double hamInitial = graph->getHam();
            double hamDiff = simPartial(*graph, xVals, yVals); //gets the change in energy if edge between nodeA and nodeB is flipped
            if(hamDiff <= 0) {  //Accepts change if it decreases energy.
                graph->flipEdge(xVals,yVals);
                graph->acceptPartial(hamDiff);
                swaps++;
            }
            else { //If swap increases energy, the change is accepted given a particular probability.
                
                double prob = exp(TINV*(-hamDiff));
                
                if(probDis(gen) <= prob) {
                    graph->flipEdge(xVals,yVals);
                    graph->acceptPartial(hamDiff);
                    increases++;
                    swaps++;
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
            
            if((sweeps % collectionTime) == 0) {
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
                    data[simNum][AVG_DEGREE].push_back(double(sum)/SIZE);
                }
                
                if(observe[EULER_CHAR]) {
                    data[simNum][EULER_CHAR].push_back(graph->getEulerChar());
                    
                }
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
    std::cout << selected << " selected." << std::endl;
    descriptor += ".csv";
    std::ofstream output(descriptor);
    output << graph->getSize() << std::endl;
    output << *graph;
    output.close();
    std::cout << std::endl;

        
    
    
    
}


bool getTF() {
    while(true) {
        char c = getchar();
        if (toupper(c) == 'Y') {
            std::cin.clear();
            std::cin.ignore(100, '\n');
            return true;
            
        }
        else if(toupper(c) == 'N') {
            std::cin.clear();
            std::cin.ignore(100, '\n');
            return false;
            
        }
        else {
            std::cout << "Invalid input." << std::endl;
            std::cout << "Please enter a 'y' or an 'n': ";
            std::cin.clear();
            std::cin.ignore(100, '\n');
            
        }
    }
}


                                                   
                                                


