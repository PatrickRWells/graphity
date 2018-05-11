//
//  graphingUtil.cpp
//  
//
//  Created by Patrick on 5/10/18.
//

#include "graphingUtil.hpp"
#include <string>
#include <iostream>
#include <cctype>

void drawMultiGraph(std::vector<std::vector<double>> xVals, std::vector<std::vector<double>> yVals) {
    std::string fileName;
    std::cout << "Input a filename for the graph: ";
    //std::cin >> fileName;
    std::cin.clear();
    std::cin.ignore(100, '\n');
    std::getline(std::cin, fileName);
    fileName += ".png";
    int numSeries = xVals.size();
    std::cout << numSeries << " data series detected" << std::endl;
    
    std::string gname;
    std::cout << "Enter a title for the graph: ";
    std::getline(std::cin, gname);
    
    bool valid = false;
    char legendIn;
    bool isLegend = false;
    
    while(!valid) {
        std::cout << "Would you like to draw a legend (y/n)?";
        std::string legendString;
        std::getline(std::cin, legendString);
        if(toupper(legendString[0]) == 'Y' ) {
            isLegend = true;
            valid = true;
        }
        
        else if(toupper(legendString[0]) == 'N') {
            valid = true;
        }
        
        else {
            std::cout << "Invalid option." << std::endl;
        }
    }
    
    TCanvas * c1 = new TCanvas("c1","Canvas",200,10,1200,800);
    TMultiGraph * multiGraph = new TMultiGraph();
    
    TGraph **graphs = new TGraph * [numSeries];
    for(int i = 0; i < numSeries; i++) {
        int numVals = xVals[i].size();
        double xArray[numVals];
        double yArray[numVals];
        for(int j = 0; j < numVals; j++) {
            xArray[j] = xVals[i][j];
            yArray[j] = yVals[i][j];
        }
        
        graphs[i] = new TGraph(numVals, xArray, yArray);
        std::string name = "gr";
        std::string title = "Graph ";
        char num = '0' + i;
        name += num;
        title += num;
        graphs[i]->SetName(name.c_str());
        graphs[i]->SetTitle(name.c_str());
        graphs[i]->SetLineColor(i+1);
        graphs[i]->SetLineWidth(1);
        multiGraph->Add(graphs[i]);
    }
    multiGraph->SetTitle(gname.c_str());
    TLegend * legend = new TLegend(.7, .9, .9, .7);
    legend->SetHeader("Legend", "C");
    if(isLegend) {
        for(int i = 0; i < numSeries; i++) {
            std::string lname;
            std::cout << "Input title of data series " << i << ": ";
            std::getline(std::cin, lname);
            legend->AddEntry(graphs[i], lname.c_str());
            
        }
    }
    
    multiGraph->Draw("AC");
    if(isLegend) {
        legend->Draw();
        
    }
    c1->Print(fileName.c_str());
    delete multiGraph;
    delete [] graphs;
    delete c1;
    
    
    
    
}
