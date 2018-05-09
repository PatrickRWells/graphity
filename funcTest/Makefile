CXX = g++
CXXFLAGS = -std=c++11
OBJS = depend/hGraph.o
HGRAPH = ../resource/hGraph/
MAIN = funcTest

INCLUDES = -I ../resource -I ../hamiltonian -I $(HOME)/root/include

default: funcTest

funcTest: depend/hGraph.o
	$(CXX) $(CXXFLAGS) $(INCLUDES) -o $(MAIN) $(OBJS) funcTest.cpp `root-config --glibs`

depend/hGraph.o: ../resource/hGraph/hGraph.cpp
	$ mkdir -p depend
	$(CXX) $(CXXFLAGS) $(INCLUDES) -c  ../resource/hGraph/hGraph.cpp -o depend/hGraph.o

clean:
	$(RM) *.o $(MAIN)
	$(RM) -r depend
