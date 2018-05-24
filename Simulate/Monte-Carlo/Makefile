CXX = g++
CXXFLAGS = -std=c++11
OBJS = depend/hGraph.o depend/graphingUtil.o depend/graphImager.o
HGRAPH = ../../resource/hGraph/
MAIN = MonteCarlo
ROOT = ../..

INCLUDES = -I $(ROOT)/resource -I $(ROOT)/hamiltonian -I $(HOME)/root/include -I $(ROOT)/lib/nauty

default: MonteCarlo

MonteCarlo: depend/hGraph.o depend/graphingUtil.o depend/graphImager.o ../../hamiltonian/hamiltonians.h MonteCarlo.cpp
	$(CXX) $(CXXFLAGS) $(INCLUDES) -o $(MAIN) $(OBJS) MonteCarlo.cpp `root-config --glibs` `pkg-config --cflags --libs cairo` $(ROOT)/lib/nauty/nauty.a

depend/graphImager.o: ../../resource/graphics/graphImager/graphImager.cpp
	$(CXX) $(CXXFLAGS) $(INCLUDES) -c ../../resource/graphics/graphImager/graphImager.cpp -o depend/graphImager.o `pkg-config --cflags --libs cairo`

depend/graphingUtil.o: ../../resource/graphics/graphUtil/graphingUtil.cpp
	$(CXX) $(CXXFLAGS) $(INCLUDES) -c ../../resource/graphics/graphUtil/graphingUtil.cpp -o depend/graphingUtil.o 


depend/hGraph.o: ../../resource/hGraph/hGraph.cpp 
	$ mkdir -p depend
	$(CXX) $(CXXFLAGS) $(INCLUDES) -c  ../../resource/hGraph/hGraph.cpp -o depend/hGraph.o

clean:
	$(RM) *.o $(MAIN)
	$(RM) -r depend
