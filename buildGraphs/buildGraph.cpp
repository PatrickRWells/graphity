//
//  buildGraph.cpp
//  
//
//  Created by Patrick on 10/14/17.
//
//

#include <iostream>
#include "hGraph.h"
#include <algorithm>
using namespace std;

int main() {
    const int SIZE = 3; //size of adjacency matrix. AKA the number of nodes in the graph
    
    int max = (SIZE*(SIZE - 1))/2; //the number of possible connections that have to be considered
    
    for(int i = 0; i < max; i++) {
        int fill[max];
        for(int l = 0; l < max; l++) { //maps bottom half of the adjacency matrix to an array of intergers that can be either 0 or 1
                                       //the 0th item of the array is is the first variable entry (1 row, 0th column) 1st item of the array is second variable entry (2 row, 0th column) 2nd item of the array is third variable entry (2 row, 1st column) and so forth;
            fill[l] = 0;
        }
        for(int j = 0; j <= i; j ++) { //loops over the number of possible 1s in the identity matrix
            fill[max - 1 - j] = 1;     //fills the array map with the correct number of 1s, starting with the last item in the array
        }

        
        do {
            MatrixXi adjMatrix = MatrixXi::Zero(SIZE, SIZE);
            int index = 0;
            for(int k = 1; k < SIZE; k++) {
                for(int m = 0; m < k; m++) {
                    adjMatrix(k, m) = fill[index];  //maps the item in the array to its corresponding entry in the bottom half of the adjacency matrix
                    adjMatrix(m, k) = fill[index];  //maps the item in the array to its corresponding entry in the top half of the adjacency matrix
                    index++;
                
                }
                
                
                
            }
            cout << adjMatrix << endl << endl;      //outputs adjacency matrix for testing purposes. Future versions will simply input the matrix into the appropriate graph object

            
            
        } while(next_permutation(fill, fill + max)); //modifies the array map to be the next possible permutation of the array for the number of 1s currently in array
                                                     //if there is no next possible permutation, the loop is terminated and the process restarts with a new number of 1s
        
        
        
    }
    
    
    
}

