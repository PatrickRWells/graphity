CXX = g++
CXXFLAGS = -std=c++11
OBJS = depend/hGraph.o depend/graphingUtil.o
HGRAPH = ../../resource/hGraph/
MAIN = expectationVals

INCLUDES = -I ../../resource -I ../../hamiltonian -I $(HOME)/root/include

default: expectationVals

expectationVals: depend/hGraph.o depend/graphingUtil.o ../../hamiltonian/hamiltonians.h expectationVals.cpp
	$(CXX) $(CXXFLAGS) $(INCLUDES) -o $(MAIN) $(OBJS) expectationVals.cpp `root-config --glibs`

depend/graphingUtil.o: ../graphics/graphUtil/graphingUtil.cpp
	$(CXX) $(CXXFLAGS) $(INCLUDES) -c ../graphics/graphUtil/graphingUtil.cpp -o depend/graphingUtil.o 


depend/hGraph.o: ../hGraph/hGraph.cpp 
	$ mkdir -p depend
	$(CXX) $(CXXFLAGS) $(INCLUDES) -c  ../hGraph/hGraph.cpp -o depend/hGraph.o

clean:
	$(RM) *.o $(MAIN)
	$(RM) -r depend
