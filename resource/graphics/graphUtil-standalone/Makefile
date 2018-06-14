CXX = g++
CXXFLAGS = -std=c++11
OBJS = graphingUtil.o
MAIN = graphCSV
INCLUDES = -I $(HOME)/root/include -I ../graphUtil

graphCSV: graphingUtil.o graphCSV.cpp
	$(CXX) $(CXXFLAGS) $(INCLUDES) -o $(MAIN) $(OBJS) graphCSV.cpp `root-config --glibs`

graphingUtil.o: ../graphUtil/graphingUtil.cpp
	$(CXX) $(CXXFLAGS) $(INCLUDES) -c ../graphUtil/graphingUtil.cpp -o graphingUtil.o `root-config --glibs`

