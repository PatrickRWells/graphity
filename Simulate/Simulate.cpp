//
//  Simulate.cpp
//  
//
//  Created by Patrick on 11/2/17.
//
//


#include <iostream>
#include <string>
#include "hGraph/hGraph.h"

using namespace std;


int main() {
    
    int size;
    char cSize;
    string filename;
    ifstream input;
    
    while(true) {
        cout << "Input graph data filename: "; //gets input file name
        cin >> filename;
        filename = "../" + filename;           //changes file reference to subdirectory
        input.open(filename);
        if(input.good()) {                  //makes sure file exists
            break;
        }
        else {
            cout << "Invalid file. Check that the filenaem is spelled correctly and that it is in the main folder." << endl;
        }
    }
    
    
    string line;
    getline(input, line);
    size = line[0] - '0'; //if a character represents an integer, subtracting character 0 from it converts its binary representation to the actual integer.
    
    int linesize = 2*size*size; //figures out how much of the line must be grabbed (including commas) to get the adjacency matrix
    int data[size*size];        //creates array that will temporarily hold adjacency matrix data;
    

    int numRead = 0;
    
    while(true) {
        getline(input,line); //gets the next line of the CSV file
        
        if(input.eof()) {  //checks if end of CSV file has been reached;
            break;
        }

        MatrixXi adjMatrix = MatrixXi::Zero(size, size);

        for(int j = 0; j < 2*size*size; j += 2) {
            adjMatrix(j/(2*size), (j/2) % size) = line[j] - '0';   //adds the integer to the adjacency matrix in the appropriate position
        }
        
        numRead++;
        
        cout << adjMatrix << endl << endl; //outputs adjacency matrix for debuggin purposes.
        
        
    }
    
    cout << "Number of graphs read " << numRead << endl;
    
    
}



