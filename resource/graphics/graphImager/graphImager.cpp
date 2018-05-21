//
//  graphImager.cpp
//  Creates an image of a graph for your viewing pleasure. Uses the cairo graphing library.
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
    /*
    The next several lines are effectively a monte-carlo simulation for the locations of the nodes on the graph.
     It minimizes a function that is based on the relative locations of the vertices.
     The energy between two connected nodes is minimized when they are a given preffered distance away from each other.
     The energy between two non-connected nodes is minimized when they are as far apart as possible.
     See Section IV in arXiv:1210.3372 [gr-qc] for details.
    
    */
    double locations[graphSize][2];
    
    for(int i = 0; i < graphSize; i++) { //generates an initialy random location for all the nodes.
        std::random_device rd; //c++11 random number generator
        std::mt19937 gen(rd());
        std::uniform_real_distribution<double> randAngle(0, 2*M_PI);
        std::uniform_real_distribution<double> randDistance(0.0, 1.0);
        //Locations are generated in polar coordinates, giving the graph a nice-looking circular shape.
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
    
    
    //Regular monte carlo simulation.
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
    
    //This is the part that actually makes the image. It uses the cairo library, which can be found at https://cairographics.org/
    cairo_surface_t *surface;
    surface = cairo_image_surface_create(CAIRO_FORMAT_ARGB32, imageSize, imageSize);
    cairo_t *cr;
    cr = cairo_create(surface);
    
    
    cairo_set_source_rgb(cr, 0, 0, 0);
    cairo_translate(cr, imageSize/2, imageSize/2 ); //sets the reference point to the center of the image.

    
    
    /*
     The generated node locations are all in a sphere of radius 1. The scale function ensures these coordinates
     are translated the correct values for the given image size, with some room for padding.
     For example, if a node was was located at x = 1, y = 0, the node would be placed 47% of the total image width
     to the right of the center of the image. The other 3% is left for the padding.
     */
    cairo_scale(cr, imageSize*0.47, imageSize*0.47);
    
    double lineWidth = 1;
    double nodeSize = 3;
    cairo_device_to_user_distance(cr, &lineWidth, &nodeSize);
    cairo_set_line_width(cr, lineWidth);
    //The scale function also scales parameters such as line width, so this ensures they are set
    //to the correct value regardless of image size.
    
    
    for(int i = 0; i < graphSize; i++) {
        for(int j = i+1; j< graphSize; j++) {
            //Lines are graphed first. The commands are fairly self explanatory.
            if(graph.isConnected(i,j)) {
                cairo_move_to(cr, locations[i][0], locations[i][1]);
                cairo_line_to(cr, locations[j][0], locations[j][1]);
                
                
            }
            
        }
    }
    //Actually draws the lines on the canvas.
    cairo_stroke(cr);
    
    //I wanted it to look pretty so set the nodes to be a nice shade of blue. Fight me.
    double r = 19;
    double g = 26;
    double b = 196;
    
    cairo_set_source_rgb(cr, r/255, g/255, b/255); //Cairo uses decimal RGB values.
    for(int i = 0; i < graphSize; i++) {
        //Cairo doesn't have an actual circle object, but it can draw arks and fill enclosed spaces.
        //This goes to the location of each individual node and drawas a 360 degree arc, then fills it.
        //Also I hate that I just said "degrees" but I can't write the pi character in Xcode so it'll have to do.
        cairo_arc(cr, locations[i][0], locations[i][1], nodeSize,  0, 2*M_PI);
        cairo_stroke_preserve(cr);
        cairo_fill(cr);

    }

     

    //Outputs the resultant figure to the file specificed by the user.
    cairo_surface_write_to_png(surface, imageName.c_str());
    std::cout << "Graph output to file " << imageName << std::endl;
    
    
    
}

double energyVal(double position[][2], hGraph graph) {
    //This function determines the "energy" for a given configuration of the nodes.
    //See Section IV in arXiv:1210.3372 [gr-qc] for details.

    int size = graph.getSize();
    double sumA = 0;
    double sumB = 0;
    double a1 = 1;
    double a2 = 1; //The prefered edge length in the unit circle based on this function turns out to be (a2/a1)^1/4.
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
    return (x2 - x1)*(x2 - x1) + (y2 - y1)*(y2 - y1); //Actually don't think this really saves any programing space but who cares. 
}

