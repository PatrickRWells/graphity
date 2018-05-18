//
//  graphImager.h
//  
//
//  Created by Patrick on 5/18/18.
//

#ifndef graphImager_h
#define graphImager_h

#include <iostream>
#include <cmath>
#include <random>
#include <cairo.h>
#include "hGraph/hGraph.h"

void graphImage(hGraph graph);
double energyVal(double position[][2], hGraph graph);
double distanceSquare(double x1, double y1, double x2, double y2);



#endif /* graphImager_h */
