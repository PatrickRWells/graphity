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


double TINV;
int SIZE;
int maxSweeps;

std::function<void(hGraph&)> simFunction;
std::function<double(hGraph&,int, int)> simPartial;




int main() {
    
    simFunction = basicSquareHam;
    simPartial = basicSquarePartial;

    
    std::vector<double> * energyStore = new std::vector<double>;
    std::vector<double> * energyStore2 = new std::vector<double>;
    
    std::cout << "How large of a graph would you like to use? ";
    std::cin >> SIZE;
    std::cout << std::endl << "Input the inverse temperature: ";
    std::cin >> TINV;
    std::cout << "How many sweeps would you like to perform? ";
    std::cin >> maxSweeps;
    
    
    hGraph * graph = new hGraph(SIZE);
    (*graph) = randomGraph(SIZE);
    
    hGraph * graphK = new hGraph(SIZE);
    *graphK = compGraph(SIZE);
    
    std::thread threadA(monteCarlo, graph, energyStore);
    monteCarlo(graphK, energyStore2);
    
    threadA.join();
    
    Int_t points = energyStore->size();
    Int_t points2 = energyStore2->size();

    Double_t x[points], y[points], x2[points2], y2[points2];


    for(int i = 0; i < points; i++) {
        x[i] = i;
        y[i] = (*energyStore)[i];
    }
    
    for(int i = 0; i < points2; i++) {
        x2[i] = i;
        y2[i] = (*energyStore2)[i];
    }
    
    TCanvas *c1 = new TCanvas("c1","Graph Examples",200,10,600,400);

    TMultiGraph *graphMulti = new TMultiGraph();
    graphMulti->SetTitle("Energy");
    
    TGraph *gr1 = new TGraph (points, x, y);
    gr1->SetLineColor(2);
    gr1->SetLineWidth(1);
    gr1->SetFillStyle(3005);
    
    TGraph *gr2 = new TGraph (points, x2, y2);
    gr2->SetLineColor(4);
    gr2->SetLineWidth(1);
    gr2->SetFillStyle(3005);

    graphMulti->Add(gr1);
    graphMulti->Add(gr2);
    
    graphMulti->Draw("AC");
    c1->Print("test.png");
    delete energyStore;
    delete energyStore2;
    
    
}

void monteCarlo (hGraph * graph, std::vector<double> * energy ) {
    std::random_device rd; //c++11 random number generator
    std::mt19937 gen(rd());
    std::uniform_int_distribution<int> randNode(0, SIZE - 1);
    std::uniform_real_distribution<> probDis(0.0, 1.0);
    
    bool run = true;
    int nodeA = 0;
    int nodeB = 0;
    int sweeps = 0;
    int sweepNUM = ((SIZE)*(SIZE-1))/2;
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
        if(n == sweepNUM) {
            n = 0;
            sweeps++;
            energy->push_back(graph->getHam());
            if((maxSweeps >= 10) && sweeps % (maxSweeps/10) == 0) {
                std::cout << sweeps << " sweeps performed" << std::endl;
            }
            
        }
        
        if(sweeps == maxSweeps) {
            run = false;
        }
        
        
        
    }
    std::cout << std::endl;
    std::cout << *graph << std::endl;
    std::cout << "Moves that increased energy: " << increases << std::endl;
    std::cout << "total swaps accepted: " << swaps <<std::endl;

    
    
    
}
