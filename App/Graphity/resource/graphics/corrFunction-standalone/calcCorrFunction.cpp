//
//  calcCorrFunction.cpp
//  
//
//  Created by Patrick on 6/14/18.
//

#include <iostream>
#include <fstream>
#include "graphingUtil.hpp"

int main() {
    std::string filename;
    std::string ofilename;
    std::ifstream input;
    std::ofstream output;
    
    while(true) {
        std::cout << "Enter the filename that contains the data for which you would like to calculate the autocorrelation function: ";
        std::getline(std::cin, filename);
        input.open(filename);
        if(input.good()) {
            break;
        }
        else {
            std::cout << "File not found. Ensure the filename and location are correct." << std::endl ;
        }
        
    }
    
    std::cout << "Enter a name for the output file: ";
    std::getline(std::cin, ofilename);
    ofilename += ".csv";
    
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
        
        data[i] = new std::vector<double> [2];;
        
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
    
    for(int i = 0; i < numLines; i++) {
        std::cout << "Calculating correlation function for data series " << i+1 << std::endl;
        correlationFn(data, i, 0, 1);
    }
    
    output.open(ofilename);
    for(int i = 0; i < numLines; i++) {
        for (int j = 0; j < data[i][1].size() - 1; j++){
            output << data[i][1][j];
            output << ",";
        }
        output << data[i][1][data[i][1].size() - 1];
        output << std::endl;
    }
    output.close();
    std::cout << "Correlation function written to file " << ofilename << std::endl;

    for(int i = 0; i < numLines; i++) {
        delete [] data[i];
    }
    delete [] data;
}
