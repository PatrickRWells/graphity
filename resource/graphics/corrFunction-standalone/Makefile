CXX = g++
CXXFLAGS = -std=c++11
OBJS = graphingUtil.o
MAIN = calcCorrFunction
INCLUDES = -I $(HOME)/root/include -I ../graphUtil

calcCorrFunction: graphingUtil.o calcCorrFunction.cpp
	$(CXX) $(CXXFLAGS) $(INCLUDES) -o $(MAIN) $(OBJS) calcCorrFunction.cpp `root-config --glibs`

graphingUtil.o: ../graphUtil/graphingUtil.cpp
	$(CXX) $(CXXFLAGS) $(INCLUDES) -c ../graphUtil/graphingUtil.cpp -o graphingUtil.o 

clean:
	$(RM) *.o $(MAIN)
