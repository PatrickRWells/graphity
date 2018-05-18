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

void drawMultiGraph(std::vector<double> data[][5], int numSeries, int observable);
void correlationFn(std::vector<double> data[][5], int run, int inData, int outData);


#ifndef graphingUtil_hpp
#define graphingUtil_hpp


#endif /* graphingUtil_hpp */
