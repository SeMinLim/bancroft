LIB = -lrt

all: $(wildcard *.cpp)
	mkdir -p bin
	g++ -o bin/bancroft_compressor $(wildcard *.cpp) -Wall -pedantic -g -O2 -lm -ffast-math
