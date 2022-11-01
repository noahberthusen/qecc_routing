CXX = g++
CXXFLAGS = -std=c++17 -fopenmp -Wall -Wextra -O3

main: main.o generator.o mat.o grid.o result.o
	$(CXX) $(CXXFLAGS) -o main main.o generator.o mat.o grid.o result.o

main.o: main.cpp mat.h generator.h result.h
	$(CXX) $(CXXFLAGS) -c main.cpp

grid.o: grid.cpp grid.h generator.h mat.h
	$(CXX) $(CXXFLAGS) -c grid.cpp

generator.o: generator.cpp generator.h
	$(CXX) $(CXXFLAGS) -c generator.cpp

result.o: result.cpp result.h
	$(CXX) $(CXXFLAGS) -c result.cpp

mat.o: mat.h
	$(CXX) $(CXXFLAGS) -c mat.cpp

clean:
	rm -rf *.o