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

void drawMultiGraph(std::vector<std::vector<double>> xVals, std::vector<std::vector<double>> yVals);
void correlationFn(std::vector<double> * data, std::vector <double> * output); 


#ifndef graphingUtil_hpp
#define graphingUtil_hpp


#endif /* graphingUtil_hpp */
