//
//  graphCSV.cpp
//  
//
//  Created by Patrick on 5/16/18.
//

#include <iostream>
#include <fstream>
#include "graphingUtil.hpp"

int main() {
    std::string filename;
    std::ifstream input;

    while(true) {
        std::cout << "Enter the filename that contains the data you would like to plot: ";
        std::getline(std::cin, filename);
        input.open(filename);
        if(input.good()) {
            break;
        }
        else {
            std::cout << "File not found. Ensure the filename and location are correct." << std::endl ;
        }
        
    }
    
    std::string line;
    int numLines;
    while(true) {
        getline(input,line); //gets the next line of the CSV file
        
        if(input.eof()) {  //checks if end of CSV file has been reached;
            break;
        }
        
        numLines++;
    }
    input.clear();
    input.seekg(0, std::ios::beg);
    std::vector<std::vector<double>> yVals;
    for(int i = 0; i < numLines; i++) {
        std::vector<double> temp;
        getline(input, line);
        int index = 0;
        char charIn;
        std::string tempStr;
        
        for(int j = 0; j < line.length(); j++) {

            charIn = line[j];
            if(charIn != ',') {
                tempStr += charIn;
            }
            else {
                temp.push_back(std::stod(tempStr));
                tempStr = "";
            }
            

        }
        temp.push_back(std::stod(tempStr));
        tempStr = "";
        yVals.push_back(temp);
        temp.clear();
        
    }
    
    input.close();
    
    std::vector<std::vector<double>>  corrFunctions;
    std::cout << "Would you like to graph the correlation functions for the given data? (y/n) ";
    std::string inString;
    std::getline(std::cin, inString);
    bool correlation = false;
    if(toupper(inString[0]) == 'Y' ) {
        std::vector<double> temp;
        correlation = true;
        std::cout << "Calculating correlation functions... " << std::endl;
        for(int i = 0; i < numLines; i++) {
            correlationFn(&yVals[i], &temp);
            corrFunctions.push_back(temp);
            temp.clear();
        }
        
    }
    
    std::vector<std::vector<double>> xVals;
    std::vector<double> temp;
    for(int i = 0; i < yVals[0].size(); i++) {
        temp.push_back(i);
    }
    for(int i = 0; i < numLines; i++) {
        xVals.push_back(temp);
    }
    
    std::cout << "Graphing data..." << std::endl;
    drawMultiGraph(xVals, yVals);
    
    if(correlation) {
        
        std::cout << "Graph correlation function..." << std::endl;
        drawMultiGraph(xVals, corrFunctions);
        
        
    }
        

    
    
}
    
    
    
    
    
    
    
    
    
    
    

