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
    
    std::vector <double> ** data = new std::vector<double> * [numLines];
    for(int i = 0; i < numLines; i++) {
        
        data[i] = new std::vector<double> [1];;
        
    }
    
    for(int i = 0; i < numLines; i++) {
        std::getline(input, line);
        line += '\n';
        int j = 0;
        std::string tempStr = "";
        while(true) {
            char ch = line[j];
            if(ch == ',') {
                data[i][0].push_back(std::stod(tempStr));
                tempStr = "";
                j++;
                continue;
            }
            else if (ch == '\n') {
                data[i][0].push_back(std::stod(tempStr));
                tempStr = "";
                j = 0;
                break;
    
            }
            else {
                tempStr += ch;
                j++;
            }
        }
    }

    std::string description[numLines];
    for(int i = 0; i < numLines; i++) {
        std::cout << "Enter a descriptor for data series " << i+1 << ": ";
        std::getline(std::cin, description[i]);
        
    }
    
    drawMultiGraph(data, numLines, 0, description, "");
    
    for(int i = 0; i < numLines; i++) {
        delete [] data[i];
    }
    delete [] data;

    
    
}
    
    
    
    
    
    
    
    
    
    
    

