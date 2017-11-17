CXX = g++
CXXFLAGS = -std=c++11
OBJS = depend/hGraph.o
HGRAPH = ../resource/hGraph/
MAIN = Simulate

INCLUDES = -I ../resource -I ../hamiltonian

default: Simulate

Simulate: depend/hGraph.o ../hamiltonian/hamiltonians.h
	$(CXX) $(CXXFLAGS) $(INCLUDES) -o $(MAIN) $(OBJS) Simulate.cpp

depend/hGraph.o: ../resource/hGraph/hGraph.cpp 
	$ mkdir -p depend
	$(CXX) $(CXXFLAGS) $(INCLUDES) -c  ../resource/hGraph/hGraph.cpp -o depend/hGraph.o

clean:
	$(RM) *.o $(MAIN)
	$(RM) -r depend
