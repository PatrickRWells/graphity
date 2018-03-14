CXX = g++
CXXFLAGS = -std=c++11
OBJS = depend/hGraph.o
HGRAPH = ../../resource/hGraph/
MAIN = stats

INCLUDES = -I ../../resource -I ../../hamiltonian

default: stats

stats: depend/hGraph.o
	$(CXX) $(CXXFLAGS) $(INCLUDES) -o $(MAIN) $(OBJS) stats.cpp

depend/hGraph.o: ../../resource/hGraph/hGraph.cpp
	$ mkdir -p depend
	$(CXX) $(CXXFLAGS) $(INCLUDES) -c  ../../resource/hGraph/hGraph.cpp -o depend/hGraph.o

clean:
	$(RM) *.o $(MAIN)
	$(RM) -r depend