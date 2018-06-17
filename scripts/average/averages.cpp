//
//  averages.cpp
//  
//
//  Created by Patrick on 5/25/18.
//

#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>

bool getTF();


int main() {
    
    
    std::string filename;
    std::ifstream input;
    
    while(true) {
        std::cout << "Enter the filename that contains the data you would like to calculate values on: ";
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
    int numLines = 1;
    getline(input,line); //gets the next line of the CSV file
    int lineLength = line.length();
    while(true) {
        getline(input,line); //gets the next line of the CSV file
        
        if(input.eof()) {  //checks if end of CSV file has been reached;
            break;
        }
        
        numLines++;
    }
    
    
    input.clear();
    input.seekg(0, std::ios::beg);
    std::vector <double> data[numLines];
    
    for(int i = 0; i < numLines; i++) {
        getline(input,line);
        std::string temp = "";
        for(int j = 0; j <= line.length(); j++) {
            if ((line[j] == ',') || (j == (line.length()))) {
                data[i].push_back(std::stod(temp));
                temp = "";
                
            }
            else {
                temp += line[j];
            }
        }
    }
    
    double averages[numLines];
    
    std::cout << "Data read..." << std::endl;
    for(int i = 0; i < numLines; i++ ) {
        double sum = 0;
        for(int j = 0; j < data[i].size(); j++) {
            sum += data[i][j];
        }
        averages[i] = sum/(data[i].size());
        std::cout << "Data set " << i+1 << " average: " << sum/(data[i].size()) << std::endl;
    }
    
    std::cout << "Would you like to calculate the average of the square of the data? (y/n) ";
    if(getTF()) {
        for(int i = 0; i < numLines; i++ ) {
            double sum = 0;
            for(int j = 0; j < data[i].size(); j++) {
                sum += (data[i][j]*data[i][j]);
            }
            std::cout << "Data set " << i+1 << " average of squares: " << sum/(data[i].size()) << std::endl;
        }
    }
    
    std::cout << "Would you like to calculate the standard deviation of the data? (y/n) ";
    if(getTF()) {
        for(int i = 0; i < numLines; i++) {
            double sum = 0;
            for(int j = 0; j  < data[i].size(); j++) {
                double partial = data[i][j] - averages[i];
                sum += partial*partial;
                
            }
            sum = sum/(data[i].size() - 1);
            sum = sqrt(sum);
            std::cout << "Data set " << i+1 << " standard deviaton: " << sum << std::endl;
            
        }
    
    
    }
    
    
    
    
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

