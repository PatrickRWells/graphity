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



bool getTF(); //Utility function. Returns true if user types 'y' and false if user types 'n'

bool WOLFF = false;
int NUM_CORES = 4; //Number of cores in the CPU
int threadsPer = 0;
double TINV;    //Beta
int SIZE;       //Number of nodes in the graphs
int maxSweeps;  //Number of sweeps that will be performed
int collectionTime = 1;

/////////--------------It is high recommended that you read the user guide entry for the monte-carlo simulation prior to editing this file----------//////////


//------Simulation Functions------//
std::function<void(hGraph&)> simFunction;
std::function<double(hGraph&,std::vector<int>,std::vector<int>)> simPartial;
//All these variables are declared to be global so that they can be used by multiple threads.


int main() {
    
    
    
    std::vector<bool> observables (NUM_OBSERVABLES, false); // true/false. Tells the program which data to save.
    std::string folder;
    std::cout << "What folder would you like to place results in? Leave blank to place results in the current working directory. ";
    std::getline(std::cin, folder);
    if(folder.length() > 0 && folder[0] != ' ') {
        if(folder[folder.length() -1] != '/') {
            folder += '/';
        }
    }
    
    else {
        folder = "";
    }

    std::string energyOutput = folder + "energy.csv";
    
    std::ofstream engOut(energyOutput); //Energy and simulation parameters will always be outputted
    if(!engOut.is_open()) {
        std::string cmd = "mkdir " + folder;
        system(cmd.c_str());
        engOut.open(energyOutput);
        if(!engOut.is_open()) {
            std::cout << "Could not create files in the given folder... Terminating" << std::endl;
            exit(3);
            
        }
        std::cout << "The folder specified does not exist and has been automatically created" << std::endl;
        
    }
    

    
    
    std::ofstream paramOut(folder + "parameters.txt");
    std::ofstream allOut;


    std::ofstream outStreams[NUM_OBSERVABLES]; //Forward declaration of output file streams for observables. These may or may not be used depending on the data the user wishes to collect

    //Sets simulation functions. This is the section that is edited by the python script. Can also be set manually
	simFunction = testHamHam;
	simPartial = testHamPartial;
    //

    std::cout << "Would you like to input a graph from a file? (y/n) "; //Checks if the user would like to input a graph from a file (vs. using randomly initialized graphs);
    observables[USER_IN] = getTF(); //The constants corresponding to parameters (such as USER_IN) are defined in graphingUtil.hpp

    
    hGraph ** graphs; //Double pointer that contains the graphs
    int numGraphs = 0;
    
    if(observables[USER_IN]) { //reads in the graph file if the user wants to
        readGraphFile(&graphs, numGraphs);
        SIZE = graphs[0]->getSize();
        std::cout << "Graph size detected to be " << SIZE << std::endl;
        std::cout << "Would you like to manually set how often to collect data? (y/n)";
        
        if(getTF()) { //Gives the user the option to pick how frequently they would like to collect data.
            std::string temp;
            std::cout << "Input how often you would like to collect data (in sweeps): ";
            std::cin >> temp;
            collectionTime = std::stoi(temp);
        }
        else {
            std::cout << "Data will be collected every sweep." << std::endl;
        }
    }
    
    
    else {
        numGraphs = 3;
    }
    std::vector<double> ** data = new std::vector<double> * [numGraphs];
    for(int i = 0; i < numGraphs; i++) {
        data[i] = new std::vector<double> [NUM_OBSERVABLES]; //Creates a 2-d array of vectors that contain storage for all observables for all graphs
        
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
    threadsPer = NUM_CORES/numGraphs;
    std::cin.clear(); //Clears the cin buffer to ensure it all works as expected later.
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

    if(observables[AVG_DEGREE]) {
        std::cout << "Would you like to calculate the correlation function for the average node degree? (y/n) ";
        observables[AVG_DEGREE_CORR] = getTF();
    }
    
    std::cout << "Would you like to calculate the Euler Characteristic for the graph? (y/n) ";
    observables[EULER_CHAR] = getTF();
    
    
    paramOut << "Graph size: " << SIZE << std::endl;
    paramOut << "Inverse Temperature: " << TINV << std::endl;
    paramOut << "Sweeps performed: " << maxSweeps << std::endl;
    paramOut << "Source term: " << SOURCE << std::endl;
    paramOut << "Sweeps between measurements: " << collectionTime << std::endl;
    paramOut.close(); //This file contians all the simulation's parameters for future reference.
    
    
    

    std::cin.clear();
    
    
    std::string description[numGraphs];
    if(observables[USER_IN]) {
        for(int i = 0; i < numGraphs; i++){
            std::cout << "Enter a short descriptor for graph " << i << ": "; //Prompts the users for ways to describe the graph if the graphs were input from a file.
            getline(std::cin, description[i]);
        }
        //If the user is inputting graphs from file, this typically means they will be calculating the dimensionality
        //As such, the simulations are performed in series to most efficiently use the processor.
        //The monteCarlo function is implemented below main.
        for(int i = 0; i < numGraphs; i++) { //Runs the monte-carlo simulations
            monteCarlo(graphs[i], observables, data, i, true, description[i]);
            
        }
        
    }

    else { //If the user did not input the graphs from a file, two random graphs and an empty graph are used as seeds.
        graphs = new hGraph * [3];
        graphs[0] = new hGraph(SIZE);
        *graphs[0] = randomGraph(SIZE, 0);
        graphs[1] = new hGraph(SIZE);
        *graphs[1] = randomGraph(SIZE, 0);
        graphs[2] = new hGraph(SIZE);
        *graphs[2] = zeroGraph(SIZE);
        description[0] = "Random Graph 1";
        description[1] = "Random Graph 2";
        description[2] = "Empty Graph";

        if(!observables[DIMEN]) { //If dimensionality is not being measured, the simulations are run in paralell
            std::thread threadA(monteCarlo, graphs[0], observables, data, 0, true, description[0]);
            std::thread threadB(monteCarlo, graphs[1], observables, data, 1, true, description[1]);
            monteCarlo(graphs[2], observables, data, 2, true, description[2]);

            threadA.join();
            threadB.join();

            
        }
        
        else { //if dimensionality IS being measured, the simulations are run in series. This is the most efficient way to use the processor.
            monteCarlo(graphs[0], observables, data, 0, true, description[0]);
            monteCarlo(graphs[1], observables, data, 1, true, description[1]);
            monteCarlo(graphs[2], observables, data, 2, true, description[2]);

            
            
            
            
        }
    }
    
    for(int i = 0; i < numGraphs; i++) {
        description[i] = folder + description[i] + ".csv";
        
        std::ofstream output(description[i]);
        output << graphs[i]->getSize() << std::endl;
        output << *graphs[i];
        output.close();

        
        
    }

    
    allOut.open(folder + "allOut.csv");
    allOut << SIZE << std::endl;
    for(int i = 0; i < numGraphs; i++) {
        allOut << *(graphs[i]); //Outputs all the graphs to a single file. Graphs are also outputted to individual files.
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
    //Correlation functions are now calculated. Implementation of correlation can be found in the graphUtil.cpp file in resources/graphics/graphUtil/
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
    
    
    //All data sets that were collected are now outputed to CSVs
    
    if(observables[DIMEN_CORR]) { //Outputs dimensionality correlation function to a CSV file
        outStreams[DIMEN_CORR].open(folder + "dimenCorrelation.csv");
        for (int i = 0; i < numGraphs; i++) {
            for(int j = 0; j < data[i][DIMEN_CORR].size(); j++) {
                outStreams[DIMEN_CORR] << data[i][DIMEN_CORR][j];
                if(j != data[i][DIMEN_CORR].size() - 1) {
                    outStreams[DIMEN_CORR] << ',';
                }
                
            }
            outStreams[DIMEN_CORR] << '\n';
        }
        outStreams[DIMEN_CORR].close();
    }

    if(observables[DIMEN]) { //outputs dimensionality to a CSV file
        outStreams[DIMEN].open(folder + "dimensionality.csv");
        for (int i = 0; i < numGraphs; i++) {
            for(int j = 0; j < data[i][DIMEN].size(); j++) {
                outStreams[DIMEN] << data[i][DIMEN][j];
                if(j != data[i][DIMEN].size()-1) {
                    outStreams[DIMEN] << ',';
                }
            }
            outStreams[DIMEN] << '\n';
        }
        outStreams[DIMEN].close();
    }
    
    if(observables[ENERGY_CORR]) { //Outputs the energy correlation function to a CSV file
        outStreams[ENERGY_CORR].open(folder + "energyCorrelation.csv");
        for (int i = 0; i < numGraphs; i++) {
            for(int j = 0; j < data[i][ENERGY_CORR].size(); j++) {
                outStreams[ENERGY_CORR] << data[i][ENERGY_CORR][j];
                if(j != data[i][ENERGY_CORR].size()-1) {
                    outStreams[ENERGY_CORR] << ',';
                }
            }
            outStreams[ENERGY_CORR] << '\n';
        }
        outStreams[ENERGY_CORR].close();
    }
    
    if(observables[AVG_DEGREE]) { //Outputs average node degree to a CSV file
        outStreams[AVG_DEGREE].open(folder + "averageDegree.csv");
        for (int i = 0; i < numGraphs; i++) {
            for(int j = 0; j < data[i][AVG_DEGREE].size(); j++) {
                outStreams[AVG_DEGREE] << data[i][AVG_DEGREE][j];
                if(j != data[i][AVG_DEGREE].size()-1) {
                    outStreams[AVG_DEGREE] << ',';
                }
            }
            outStreams[AVG_DEGREE] << '\n';
        }
        outStreams[AVG_DEGREE].close();
    }
    
    if(observables[AVG_DEGREE_CORR]) { //Outputs average degree correlation function to CSV file
        outStreams[AVG_DEGREE_CORR].open(folder + "averageDegreeCorrelation.csv");
        for (int i = 0; i < numGraphs; i++) {
            for(int j = 0; j < data[i][AVG_DEGREE_CORR].size(); j++) {
                outStreams[AVG_DEGREE_CORR] << data[i][AVG_DEGREE_CORR][j];
                if(j != data[i][AVG_DEGREE_CORR].size()-1) {
                    outStreams[AVG_DEGREE_CORR] << ',';
                }
            }
            outStreams[AVG_DEGREE_CORR] << '\n';
        }
        outStreams[AVG_DEGREE_CORR].close();
    }

    if(observables[EULER_CHAR]) { //Outputs  euler characteristic to CSV file
        outStreams[EULER_CHAR].open(folder + "eulerChar.csv");
        for (int i = 0; i < numGraphs; i++) {
            for(int j = 0; j < data[i][EULER_CHAR].size(); j++) {
                outStreams[EULER_CHAR] << data[i][EULER_CHAR][j];
                if(j != data[i][EULER_CHAR].size()-1) {
                    outStreams[EULER_CHAR] << ',';
                }
            }
            outStreams[EULER_CHAR] << '\n';
        }
        outStreams[EULER_CHAR].close();
    }
    
    
    
    
    std::cout << "Data outputed to CSV files" << std::endl;
    
    
    std::cin.clear();
    
    std::cout << "Would you like to plot the energy? (y/n) ";
    
    if(getTF()) {
        std::cout << "Plotting energy..." << std::endl; //plots energy
        drawMultiGraph(data, numGraphs, ENERGY, description, folder); //See the graphics folder in resources for information, or the corresponding entry in the user guide
    }
    
    if(observables[ENERGY_CORR]) {
        std::cout << "Would you like to plot the energy correlation function? (y/n) ";
        if(getTF()) {
            std::cout << "Plotting energy correlation function... " << std::endl; //plots energy correlation function
            drawMultiGraph(data, numGraphs, ENERGY_CORR, description, folder);
        }
        
    }

    if(observables[DIMEN]) {
        std::cout << "Would you like to plot the dimensionality? (y/n) ";
        if(getTF() ){
            std::cout << "Plotting dimensionality... " << std::endl; //plots dimensionality.
            drawMultiGraph(data, numGraphs, DIMEN, description, folder);
        }
    }
    
    if(observables[DIMEN_CORR]) {
        std::cout << "Would you like to plot the dimensionality correlation function? (y/n) ";
        if(getTF()) {
            std::cout << "Plotting dimensionality correlation function... " << std::endl; //plots dimensionality correlation function
            drawMultiGraph(data, numGraphs, DIMEN_CORR, description, folder);
        }
    }
    
    if(observables[AVG_DEGREE]) {
        std::cout << "Would you like to plot the average node degree? (y/n) ";
        if(getTF()) {
            std::cout << "Plotting average node degree" << std::endl; //plots average node degree
            drawMultiGraph(data, numGraphs, AVG_DEGREE, description, folder);
            
        }
        
    }
    
    if(observables[AVG_DEGREE_CORR]) {
        std::cout << "Would you like to plot the average node degree correlation function? (y/n) ";
        if(getTF()) {
            std::cout << "Graph average node degree correlation function" << std::endl; //plots dimensionality correlation function
            drawMultiGraph(data, numGraphs, AVG_DEGREE_CORR, description, folder);
            
        }
        
    }
    
    if(observables[EULER_CHAR]) { //Plots euler characteristic
        std::cout << "Would you like to plot the Euler characteristic? (y/n) ";
        if(getTF()) {
            std::cout << "Graphing Euler Characteristic" << std::endl;
            drawMultiGraph(data, numGraphs, EULER_CHAR, description, folder);
        }
        
        
    }
    
    std::cout << "Would you like to make images of the graphs? (y/n) "; //Asks the user if they want to create a fancy PNG of the graph. See graphImager in resources/graphics
    if(getTF()) {
        for(int i = 0; i < numGraphs; i++) {
            std::cout << "Imaging graph " << description[i] << "..." << std::endl;
            graphImage(*graphs[i], folder);
        }
    }

    for(int i = 0; i < numGraphs; i++) {
        delete graphs[i];
    }
    delete [] graphs;

    for(int i = 0; i < numGraphs; i++) {
        delete [] data[i];
    }

                           
    delete [] data;
    std::cout << "If you plan on putting future results in the " << folder << " folder, be sure to move the data elsewhere to avoid overwriting it." << std::endl;
    
}
                           

void monteCarlo (hGraph * graph, std::vector<bool> observe, std::vector<double> ** data, int simNum, bool progress, std::string descriptor) {
    int selected = 0;
    graph->setThreads(threadsPer);
    std::random_device rd; //c++11 random number generator
    std::mt19937 gen(rd());
    std::uniform_int_distribution<int> randNode(0, SIZE - 1); //Generates a random integer corresponding to anode
    std::uniform_real_distribution<> probDis(0.0, 1.0); //generates a number between 0 and 1 for probability purposes
    std::uniform_int_distribution<int> randNum(1, ceil(double(SIZE)/2)); //Generates a number between 1 and NUM_NODES/2 for use in determining how many edges to flips
    //Values needed for individual simulation steps.
    bool run = true;
    int nodeA = 0;
    int nodeB = 0;
    int sweeps = 0;
    int sweepNUM = ((SIZE)*(SIZE-1))/2;//Number of edge flips in a single sweep;
    int n = 0;
    simFunction(*graph); //calculates inital graph energy;
    data[simNum][ENERGY].push_back(graph->getHam());
    
    if(observe[DIMEN]) { //Calculates initial dimensionality (if requested)
        data[simNum][DIMEN].push_back(graph->getDimension());
    }
    
    if(observe[AVG_DEGREE]) { //Calculates initial average node degree (if requested)
        double sum = 0;
        for(int i = 0; i < SIZE; i++) {
            sum += graph->getDegree(i);
        }
        data[simNum][AVG_DEGREE].push_back(sum/SIZE);
    }

    
    while(run) { //Main loop of the monte carlo simulation.
        std::vector<int> xVals;
        std::vector<int> yVals;
        
        int numFlips = randNum(gen); //Determines how many edges will be flipped
        int numAccept = 0;
        
        bool isFull = false;
        while(!isFull) { //Generates random nodes until the number of edges determined above has been reached
            while(nodeA == nodeB) { //Randomly selects two nodes;
                nodeA = randNode(gen);
                nodeB = randNode(gen);
            }
            bool same = false;
            for(int l = 0; l < xVals.size(); l++) {
                if(nodeA == xVals[l]) { //Checks to ensure that the edge generated in the previous step has not already been picked
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
            
            if(!same) { //Adds the edg to the vectors
                xVals.push_back(nodeA);
                yVals.push_back(nodeB);
                numAccept++;
            }
            if(numAccept == numFlips) { //Terminates the loop when the number of edges to be flipped has been reached.
                isFull = true;
            }
            nodeA = 0;
            nodeB = 0;
        }
        bool select = true;
        hGraph temp(graph->getSize());
        temp = *graph;
        temp.flipEdge(xVals, yVals);
        //Creates a second hGraph. Equals the original graph with the edges generated above flipped.


        mut.lock(); //Because of an implementation choice in Nauty/Traces, only one thread can use it at a time. Mut is a global mutex that ensures this is the case.
        bool iso = isIsomorphic(*graph, temp);
        mut.unlock();

        
        if (iso) {
            continue;         //If the two graphs are isomorphic, the loop entirely starts over.

        }
        
        //If they are not isomorphic, Calculate the automorphism group sizes. If the size of new is greater than size of old, the new graph can be tested
        //If the size of old is greater than size of the new, generate a random number. Test the new graph with probability (new group size)/(old group size)

    
        else {
            mut.lock();
            std::vector<double> aGrp1 = graph->autoGroupSize();
            std::vector<double> aGrp2 = temp.autoGroupSize();
            mut.unlock();
            double randAuto = (aGrp2[0]/aGrp1[0])*pow(10, aGrp2[1] - aGrp1[1]);

            
            if(randAuto >= 1.0) {
                selected++;
                select = true;
            }
            else {
                
                if(probDis(gen) <= randAuto) {
                    selected++;
                    select = true;
                }
            }
            
        }
        

        if(select) {
            double hamInitial = graph->getHam();
            double hamDiff = simPartial(*graph, xVals, yVals); //gets the change in energy if edges between nodeA and nodeB is flipped
            if(hamDiff <= 0) {  //Accepts change if it decreases energy.
                graph->flipEdge(xVals,yVals);
                graph->acceptPartial(hamDiff);
            }
            else { //If swap increases energy, the change is accepted given a particular probability.
                
                double prob = exp(TINV*(-hamDiff));
                
                if(probDis(gen) <= prob) {
                    graph->flipEdge(xVals,yVals);
                    graph->acceptPartial(hamDiff);
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
                
                           
                if(observe[DIMEN]) { //Calculates and stores dimension if requested
                    graph->setThreads(threadsPer);
                    data[simNum][DIMEN].push_back(graph->getDimension());
                    
                }
                
                if(observe[AVG_DEGREE]) { //Calculates and stores average degree if requested
                    data[simNum][AVG_DEGREE].push_back(graph->getAvgDegree());
                }
                
                if(observe[EULER_CHAR]) { //Calculats and stores Euler Characteristic if requested
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


void wolffAlgorithm(hGraph * graph, std::vector<bool> observe, std::vector<double> data[][NUM_OBSERVABLES], int simNum, bool progress, std::string descriptor) {
    std::cout << *graph << std::endl;
    int selected = 0;
    graph->setThreads(threadsPer);
    std::random_device rd; //c++11 random number generator
    std::mt19937 gen(rd());
    std::uniform_int_distribution<int> randNode(0, SIZE - 1);
    std::uniform_real_distribution<> probDis(0.0, 1.0);
    
    bool run = true;
    int sweeps = 0;
    int sweepNUM = ((SIZE)*(SIZE-1))/2;//Number of edge flips in a single sweep;
    int n = 0;
    int swaps = 0;      //number of times an edge flip is accepted.
    
    simFunction(*graph); //calculates inital graph energy;
    data[simNum][ENERGY].push_back(graph->getHam());
    
    if(observe[DIMEN]) {
        data[simNum][DIMEN].push_back(graph->getDimension());
    }
    
    if(observe[AVG_DEGREE]) {
        data[simNum][AVG_DEGREE].push_back(graph->getAvgDegree());
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
                                                


