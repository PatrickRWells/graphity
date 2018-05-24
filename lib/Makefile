CXX = g++
CXXFLAGS = -std=c++11
OBJS = depend/hGraph.o nauty/nauty.a
HGRAPH = ../resource/hGraph/
MAIN = nautyTest

INCLUDES = -I ../resource -I ../hamiltonian -I nauty

default: nautyTest

nautyTest: depend/hGraph.o
	$(CXX) $(CXXFLAGS) $(INCLUDES) -o $(MAIN) $(OBJS) nautyTest.cpp

depend/hGraph.o:
	$ mkdir -p depend
	$(CXX) $(CXXFLAGS) $(INCLUDES) -c ../resource/hGraph/hGraph.cpp -o depend/hGraph.o

clean:
	$(RM) *.o $(MAIN)
	$(RM) -r depend