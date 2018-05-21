//
//  graphingUtil.hpp
//  
//
//  Created by Patrick on 5/10/18.
//

#include "TGraphErrors.h"
#include "TCanvas.h"
#include "TMultiGraph.h"
#include "TLegend.h"
#include "TAxis.h"

const int NUM_OBSERVABLES = 6;
/*
 This is a very special constant. This is the only file that is referenced by all the
 individual pieces responsible for collecting, managing, and plotting Monte-Carlo Data
 This number tells the compiler how big to make some arrays. Just change it if you wish to add more
 And then add onto the list of constants below. 
 */

const int DIMEN = 0;
const int DIMEN_CORR = 1; //These constants determine which vector from the array
const int ENERGY_CORR = 2; //declared below to look at when looking for a particular data set
const int USER_IN = 3;
const int ENERGY = 4;
const int AVG_DEGREE = 5;


void drawMultiGraph(std::vector<double> data[][NUM_OBSERVABLES], int numSeries, int observable);
void correlationFn(std::vector<double> data[][NUM_OBSERVABLES], int run, int inData, int outData);


#ifndef graphingUtil_hpp
#define graphingUtil_hpp


#endif /* graphingUtil_hpp */
