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
    std::cout << "Input a filename for the plot: ";
    //std::cin >> fileName;
    std::getline(std::cin, fileName);
    fileName += ".png";
    int numSeries = xVals.size();
    int values = xVals[0].size();
    
    std::string lowBoundS;
    std::string highBoundS;
    
    int lowBound;
    int highBound;
    
    std::cout << numSeries << " data series detected" << std::endl;
    std::cout << "After how many sweeps would you like to start plotting? [0.." << values-1 << "] ";
    std::getline(std::cin, lowBoundS);
    lowBound = std::stoi(lowBoundS);
    
    std::cout << "After how many sweeps would you like to stop plotting? [" << lowBound+1 << ".." << values-1 << "] ";
    std::getline(std::cin, highBoundS);
    highBound = std::stoi(highBoundS);


    std::string gname;
    std::cout << "Enter a title for the plot: ";
    std::getline(std::cin, gname);
    
    std::string xname;
    std::cout << "Enter a title for the x-axis: ";
    std::getline(std::cin, xname);
    
    std::string yname;
    std::cout << "Enter a title for the y-axis: ";
    std::getline(std::cin, yname);
    
    bool valid = false;
    char legendIn;
    bool isLegend = false;
    
    while(!valid) {
        std::cout << "Would you like to draw a legend (y/n)? ";
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
        int numValues = highBound - lowBound + 1;
        double xArray[numValues];
        double yArray[numValues];
        for(int j = lowBound; j <= highBound; j++) {
            
            xArray[j-lowBound] = xVals[i][j];
            yArray[j-lowBound] = yVals[i][j];
            
            
        }
        
        
        graphs[i] = new TGraph(numValues, xArray, yArray);
        std::string name = "gr";
        std::string title = "Graph ";
        char num = '0' + i;
        name += num;
        title += num;
        graphs[i]->SetName(name.c_str());
        graphs[i]->SetTitle(name.c_str());
        graphs[i]->SetMarkerColor(i+1);
        graphs[i]->SetMarkerStyle(8);
        graphs[i]->SetMarkerSize(0.5);
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
    
    multiGraph->Draw("AP");
    if(isLegend) {
        legend->Draw();
        
    }
    multiGraph->GetXaxis()->SetTitle(xname.c_str());
    multiGraph->GetYaxis()->SetTitle(yname.c_str());

    c1->Print(fileName.c_str());
    delete multiGraph;
    delete [] graphs;
    delete c1;
    
    
    
    
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

