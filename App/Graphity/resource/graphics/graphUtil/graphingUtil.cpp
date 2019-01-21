//
//  graphingUtil.cpp
//  It's actually a plotting tool, but I'm too lazy to change the name and then edit makefiles.
//
//  Created by Patrick on 5/10/18.
//

#include "graphingUtil.hpp"
#include <string>
#include <iostream>
#include <cctype>
#include <array>
#include <cstdarg>




void drawMultiGraph(std::vector<double> ** data, int numSeries, int observable, std::string descriptors[], std::string folder) {
    /*
     Takes in three parameters. A double array with all the data. An integer indicating the number
     of data series (typically the number of graphs simulated over in the Monte-Carlo simulation,
     and one of the integers defined in graphingUtil.hpp, which tells it which observable you are actually
     graphing.
     */
    

    if(folder.size() != 0 && folder[0] != ' ') {
        if(folder[folder.size() -1] != '/') {
            folder += '/';
        }
    }
    else {
        folder = "";
    }
    
    std::string fileName;
    std::cout << "Input a filename for the plot: ";
    std::getline(std::cin, fileName);
    fileName += ".png";
    fileName = folder + fileName;
    int values = data[0][observable].size();
    
    std::string lowBoundS;
    std::string highBoundS;
    
    int lowBound;
    int highBound;
    //Allows the user to decide over what range of the values they would like to plot
    std::cout << numSeries << " data series detected" << std::endl;
    std::cout << "After how many sweeps would you like to start plotting? [0.." << values-1 << "] ";
    std::getline(std::cin, lowBoundS);
    lowBound = std::stoi(lowBoundS);
    
    std::cout << "After how many sweeps would you like to stop plotting? [" << lowBound+1 << ".." << values-1 << "] ";
    std::getline(std::cin, highBoundS);
    highBound = std::stoi(highBoundS);


    //These are self explanatory.
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
    
    //Asks if the user wants to draw a legend with the plot.
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
    
    //Uses Cern Root for the actual plotting.
    TCanvas * c1 = new TCanvas("c1","Canvas",200,10,1200,800);
    TMultiGraph * multiGraph = new TMultiGraph();
    
    //A "TGraph" is created for each plot, which are then all combined into a "TMultiGraph"
    TGraph **graphs = new TGraph * [numSeries];
    for(int i = 0; i < numSeries; i++) {
        int numValues = highBound - lowBound + 1;
        double xArray[numValues];
        double yArray[numValues];
        for(int j = lowBound; j <= highBound; j++) {
            //The TGraph object takes two arrays as its input.
            xArray[j-lowBound] = j;
            yArray[j-lowBound] = data[i][observable][j];
            //If the low bound is not the first data point, things are shifted.
            
        }
        
        //Creates the "TGraph" object with the data points.
        graphs[i] = new TGraph(numValues, xArray, yArray);
        std::string name = "gr";
        std::string title = "Graph ";
        char num = '0' + i;
        name += num;
        title += num;
        graphs[i]->SetName(name.c_str());
        graphs[i]->SetTitle(name.c_str());
        //Hopefully we never run out of colors.
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
            legend->AddEntry(graphs[i], descriptors[i].c_str());
            
        }
    }
    
    multiGraph->Draw("AP"); //"AP" means the axis are drawn and the points are not conencted by a line.
    if(isLegend) {
        legend->Draw();
        
    }
    multiGraph->GetXaxis()->SetTitle(xname.c_str());
    multiGraph->GetYaxis()->SetTitle(yname.c_str());

    c1->Print(fileName.c_str()); //Saves it to file. 
    delete multiGraph;
    delete [] graphs;
    delete c1;
    
    
    
    
}

void correlationFn(std::vector<double> **data, int run, int inData, int outData) {
    //Calculates the autocorrelation function when passed a vector with data and a vector to place results
    //(pointers). See Monte Carlo Methods in Statistical Physics - 3.21
    //The variables run and inData correspond to the indices of the double array "data" that the vector with the data.
    //Run + outData corresponds to the vector that the data will be stored in.
    int tMax = data[run][inData].size();
    for(int t = 0; t < tMax; t++ ) {    //The algorithm uses the same variable names (with tp for t') as the algorithm in the book does.
        double sum = 0;                 //Calculates the function at all times.
        double partialSum = 0;
        double partialSum2 = 0;
        for(int tp = 0; tp < tMax - t; tp++) {
            partialSum += (data[run][inData][tp])*(data[run][inData][tp + t]);
        }
        
        partialSum /= (tMax - t);
        sum += partialSum;
        partialSum = 0;
        
        for(int tp = 0; tp < tMax - t; tp++) {
            partialSum += data[run][inData][tp];
        }
        partialSum /= (tMax - t);
        
        for(int tp = 0; tp < tMax - t; tp++) {
            partialSum2 += data[run][inData][tp +t];
        }
        partialSum2 /= (tMax -t);
        
        partialSum *= partialSum2;
        sum -= partialSum;
        if(t > 0) {
            sum /= data[run][outData][0];
        }
        data[run][outData].push_back(sum);
        
    }
    data[run][outData][0] = 1;
}

