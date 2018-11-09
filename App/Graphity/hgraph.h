#ifndef HGRAPH_H
#define HGRAPH_H

#include <iostream>
#include <fstream>
#include <vector>
#include <random>
#include <thread>
#include <future>
#include <limits>
#include <cmath>

#include "eigen/Eigen/Dense"
#include "lib/nauty26r11/nauty.h"

using Eigen::MatrixXi;

class absHamiltonian;

class hGraph {

        private:

            int NUM_NODES;
            MatrixXi _adjMatrix;
            int _numThreads = 1;
            double _hamiltonian;
            std::vector<int> _degVector;
            MatrixXi _paths;



            int _eulerChar = 0;
            double _dimension;
            std::vector<double> _spectralDimen;
            std::vector<double> _hausdorffDimen;
            std::vector<int> _eccentricities;
            std::vector <int> _numCliques;
            std::vector <double> curvatureAtNode;


            bool cliquesFound = false;
            bool dimensionFound = false;


            void toStream(std::ostream &os) const;
            void toFile(std::ofstream &fs) const;
            void calcDimension();
            void oldCalcDimension(); //Multithreaded version
            double oldDimension(int a, int b, bool multi);     //Single threaded version for recursive calls
            double dimension(MatrixXi amat, double *** data, std::vector<int> * trp, int lowerBound, int upperBound, bool init);
            void calculateEccen();
            void calculateHausDimen();

            void calcEulerChar();
            void calcSpectralDimen();

            void countCliques();
            void cliqueSearch(std::vector<int> R, std::vector<int> P);




        public:
            hGraph(int size, MatrixXi adjMatrix);
            hGraph(int size);
            hGraph();
            ~hGraph();

            void setMatrix(int size, MatrixXi data);
            void setHamiltonian(double val);
            void accept(absHamiltonian &ham);
            void acceptPartial(double partial);
            void setThreads(int threads);
            hGraph compliment() const;
            hGraph unitSphere(int node) const;
            void flipEdge(int nodeA, int nodeB);
            void flipEdge(std::vector<int> nodeA, std::vector<int> nodeB);
            bool equals(hGraph &rhs);



            std::vector<int> getFractionalDimen();
            double getAvgDegree() const;
            std::vector<int> getEccentricity();
            int getDiameter();
            std::vector<double> autoGroupSize() const;
            MatrixXi getShortestPaths();
            double getDimension();
            std::vector<double> getSpectralDimen();
            std::vector<double> getHausdorffDimen();
            int getDegree(int node) const;
            int getSize() const;
            int getEulerChar();
            double curvatureAt(int node);
            double getHam() const;

            bool isConnected(int row, int column) const;
            std::vector<int> numCliques();

            bool operator == (hGraph & rhs) const;
            bool operator != (hGraph &rhs) const;
            friend std::ostream &operator << (std::ostream &os, const hGraph &rhs);
            friend std::ofstream &operator <<(std::ofstream &fs, const hGraph &rhs);

};

//Forward declaration of resource functions

void readGraphFile(hGraph *** graphs, int &num);

MatrixXi unitSphere(MatrixXi matrix, int node);
void removeColumn(Eigen::MatrixXi& matrix, unsigned int colToRemove);
void removeRow(Eigen::MatrixXi& matrix, unsigned int rowToRemove);
bool isIsomorphic(hGraph graph1, hGraph graph2);
hGraph randomGraph(int size, double fillFrac);
hGraph compGraph(int size);
hGraph zeroGraph(int size);
std::vector<int> fractionAdd(std::vector<int> fractA, std::vector<int> fractB);
void fractReduce(std::vector<int> &frac);
int gcd(int numA, int numB);

#endif // HGRAPH_H
