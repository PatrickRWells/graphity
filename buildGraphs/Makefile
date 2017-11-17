CXX = g++
CXXFLAGS = -std=c++11
OBJS = depend/hGraph.o
HGRAPH = ../resource/hGraph/
MAIN = buildGraphs

INCLUDES = -I ../resource

default: buildGraphs

buildGraphs: depend/hGraph.o
	$(CXX) $(CXXFLAGS) $(INCLUDES) -o $(MAIN) $(OBJS) buildGraph.cpp

depend/hGraph.o: ../resource/hGraph/hGraph.cpp
	$ mkdir -p depend
	$(CXX) $(CXXFLAGS) $(INCLUDES) -c  ../resource/hGraph/hGraph.cpp -o depend/hGraph.o

clean:
	$(RM) *.o $(MAIN)
	$(RM) -r depend
