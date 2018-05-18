CXX = g++
CXXFLAGS = -std=c++11
OBJS = depend/hGraph.o
MAIN = graphImager
INCLUDES = -I ../../../resource -I ../../../hamiltonian

graphImager: depend/hGraph.o graphImager.cpp
	$(CXX) $(CXXFLAGS) $(INCLUDES) -o $(MAIN) $(OBJS) graphImager.cpp `pkg-config --cflags --libs cairo`

depend/hGraph.o: ../../../resource/hGraph/hGraph.cpp 
	$ mkdir -p depend
	$(CXX) $(CXXFLAGS) $(INCLUDES) -c  ../../../resource/hGraph/hGraph.cpp -o depend/hGraph.o

clean:
	$(RM) *.o $(MAIN)
	$(RM) -r depend
