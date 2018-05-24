CXX = g++
CXXFLAGS = -std=c++11
OBJS = depend/hGraph.o depend/graphingUtil.o ../lib/nauty/nauty.a
HGRAPH = ../resource/hGraph/
MAIN = funcTest

INCLUDES = -I ../resource -I ../hamiltonian -I ../lib/nauty -I $(HOME)/root/include

default: funcTest

funcTest: depend/hGraph.o depend/graphingUtil.o funcTest.cpp
	$(CXX) $(CXXFLAGS) $(INCLUDES) -o $(MAIN) $(OBJS) funcTest.cpp `root-config --glibs`

depend/graphingUtil.o: ../resource/graphics/graphUtil/graphingUtil.cpp
	$(CXX) $(CXXFLAGS) $(INCLUDES) -c ../resource/graphics/graphUtil/graphingUtil.cpp -o depend/graphingUtil.o


depend/hGraph.o: ../resource/hGraph/hGraph.cpp
	$ mkdir -p depend
	$(CXX) $(CXXFLAGS) $(INCLUDES) -c  ../resource/hGraph/hGraph.cpp -o depend/hGraph.o


clean:
	$(RM) *.o $(MAIN)
	$(RM) -r depend
