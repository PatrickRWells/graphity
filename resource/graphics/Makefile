CXX = g++
CXXFLAGS = -std=c++11
OBJS = depend/hGraph.o
HGRAPH = ../../resource/hGraph/
MAIN = graphImage
ROOT = ../..

INCLUDES = -I $(ROOT)/resource -I $(ROOT)/hamiltonian -I $(HOME)/root/include -I $(ROOT)/lib/nauty

default: graphImage

graphImage: depend/hGraph.o graphImager.cpp
	$(CXX) $(CXXFLAGS) $(INCLUDES) -o $(MAIN) $(OBJS) graphImager.cpp `pkg-config --cflags --libs cairo` $(ROOT)/lib/nauty/nauty.a

depend/hGraph.o: ../../resource/hGraph/hGraph.cpp 
	$ mkdir -p depend
	$(CXX) $(CXXFLAGS) $(INCLUDES) -c  ../../resource/hGraph/hGraph.cpp -o depend/hGraph.o

clean:
	$(RM) *.o $(MAIN)
	$(RM) -r depend
