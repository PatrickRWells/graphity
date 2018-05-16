CXX = g++
CXXFLAGS = -std=c++11
OBJS = graphingUtil.o
MAIN = graphCSV
INCLUDES = -I $(HOME)/root/include

graphCSV: graphingUtil.o graphCSV.cpp
	$(CXX) $(CXXFLAGS) $(INCLUDES) -o $(MAIN) $(OBJS) graphCSV.cpp `root-config --glibs`

graphingUtil.o: graphingUtil.cpp
	$(CXX) $(CXXFLAGS) $(INCLUDES) -c graphingUtil.cpp -o graphingUtil.o `root-config --glibs`

