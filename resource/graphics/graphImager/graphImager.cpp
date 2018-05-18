//
//  graphImager.cpp
//  
//
//  Created by Patrick on 5/17/18.
//

#include "graphImager.h"

void graphImage(hGraph graph) {
    
    double TINV = 2.5;
    int imageSize;
    std::string imageName;
    std::cout << "Input image size (in pixels) ";
    std::cin >> imageSize;
    std::cin.clear();
    std::cin.ignore(100, '\n');
    std::cout << "Input image name: ";
    getline(std::cin, imageName);
    imageName += ".png";
    
    int numSweeps = 1000;
    int num = 0;
    std::cout << "Generating image..." << std::endl;
    int graphSize = graph.getSize();

    double locations[graphSize][2];
    
    for(int i = 0; i < graphSize; i++) {
        std::random_device rd; //c++11 random number generator
        std::mt19937 gen(rd());
        std::uniform_real_distribution<double> randAngle(0, 2*M_PI);
        std::uniform_real_distribution<double> randDistance(0.0, 1.0);
        
        double tempAngle = randAngle(gen);
        double tempDistance = randDistance(gen);
        locations[i][0] = tempDistance*cos(tempAngle);
        locations[i][1] = tempDistance*sin(tempAngle);
        }
    
    int sweepCount = 0;
    int counter = 0;
    int increases = 0;
    
    std::random_device rand;
    std::mt19937 gen(rand());
    
    std::uniform_int_distribution<int> randNode(0, graphSize -1);
    std::uniform_real_distribution<double> randAngle(0, 2*M_PI);
    std::uniform_real_distribution<double> randDistance(0.0, 1.0);
    
    while(sweepCount < numSweeps) {
        double energyOld = energyVal(locations, graph);
        int node = randNode(gen);
        double xOld = locations[node][0];
        double yOld = locations[node][1];
        double angle = randAngle(gen);
        double dist = randDistance(gen);
        
        locations[node][0] = dist*cos(angle);
        locations[node][1] = dist*sin(angle);
        double energyNew = energyVal(locations, graph);
        
        if(energyNew > energyOld) {
            double prob = exp((TINV)*(energyOld - energyNew));
            double randProb = randDistance(gen);
            if(randProb > prob) {
                locations[node][0] = xOld;
                locations[node][1] = yOld;
            }
            
            else {
                increases++;
            }
        }
        
        counter++;
        if(counter == graphSize) {
            sweepCount++;
            counter = 0;
            
        }
        
    }
    
    cairo_surface_t *surface;
    surface = cairo_image_surface_create(CAIRO_FORMAT_ARGB32, imageSize, imageSize);
    cairo_t *cr;
    cr = cairo_create(surface);
    cairo_set_source_rgb(cr, 0, 0, 0);
    cairo_translate(cr, imageSize/2, imageSize/2 );

    cairo_scale(cr, imageSize*0.45, imageSize*0.45);
    double x = 1;
    double y = 3;
    cairo_device_to_user_distance(cr, &x, &y);
    cairo_set_line_width(cr, x);
    
    
    
    for(int i = 0; i < graphSize; i++) {
        for(int j = i+1; j< graphSize; j++) {
         
            if(graph.isConnected(i,j)) {
                cairo_move_to(cr, locations[i][0], locations[i][1]);
                cairo_line_to(cr, locations[j][0], locations[j][1]);
                
                
            }
            
        }
    }
    
    cairo_stroke(cr);
    
    double r = 19;
    double g = 26;
    double b = 196;
    
    cairo_set_source_rgb(cr, r/255, g/255, b/255);
    for(int i = 0; i < graphSize; i++) {
    
        cairo_arc(cr, locations[i][0], locations[i][1], y,  0, 2*M_PI);
        cairo_stroke_preserve(cr);
        cairo_fill(cr);

    }

     

    
    cairo_surface_write_to_png(surface, imageName.c_str());
    std::cout << "Graph output to file " << imageName << std::endl;
    
    
    
}

double energyVal(double position[][2], hGraph graph) {
    int size = graph.getSize();
    double sumA = 0;
    double sumB = 0;
    double a1 = 1;
    double a2 = 1;
    double a3 = 5;
    for(int i = 0; i < size; i++) {
        for (int j = i+1; j < size; j++) {
            if(graph.isConnected(i,j)) {
                double partial = 0;
                double distanceSq = distanceSquare(position[i][0], position[i][1], position[j][0], position[j][1]);
                partial += a1*distanceSq;
                partial += a2/distanceSq;
                sumA += partial;
            }
        }
    }
    
    for(int i = 0; i < size; i++) {
        for(int j = i+1; j < size; j++) {
            if(!graph.isConnected(i, j)) {
                double distanceSq = distanceSquare(position[i][0], position[i][1], position[j][0], position[j][1]);
                sumB += a3/distanceSq;
            }
        }
    }
    
    return sumA + sumB;
    
}

double distanceSquare(double x1, double y1, double x2, double y2) {
    return (x2 - x1)*(x2 - x1) + (y2 - y1)*(y2 - y1);
}

