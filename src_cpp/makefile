CXX = g++
CXXFLAGS = -std=c++17 -fopenmp -Wall -Wextra -O3

main: main.o generator.o mat.o grid.o result.o circle.o point.o
	$(CXX) $(CXXFLAGS) -o main main.o generator.o mat.o grid.o result.o circle.o point.o

main.o: main.cpp mat.h generator.h result.h circle.h
	$(CXX) $(CXXFLAGS) -c main.cpp

grid.o: grid.cpp grid.h generator.h mat.h
	$(CXX) $(CXXFLAGS) -c grid.cpp

generator.o: generator.cpp generator.h point.h
	$(CXX) $(CXXFLAGS) -c generator.cpp

result.o: result.cpp result.h
	$(CXX) $(CXXFLAGS) -c result.cpp

circle.o: circle.cpp circle.h point.h
	$(CXX) $(CXXFLAGS) -c circle.cpp

mat.o: mat.cpp mat.h
	$(CXX) $(CXXFLAGS) -c mat.cpp

point.o: point.cpp point.h
	$(CXX) $(CXXFLAGS) -c point.cpp

clean:
	rm -rf *.o