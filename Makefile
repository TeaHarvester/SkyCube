CC = g++
CXXFLAGS = -Wall -g -O2 -fopenmp -IInclude

ECOSIM: Main.o Map.o Atmosphere.o Earth.o SkyCubeMP.o SkyCube.o
	$(CC) $(CXXFLAGS) Main.o Map.o Atmosphere.o Earth.o SkyCubeMP.o SkyCube.o -o $@

Main.o: Source/Main.cpp Source/Map.cpp 
	$(CC) $(CXXFLAGS) -c Source/Main.cpp -o Main.o

Map.o: Source/Map.cpp Source/Atmosphere.cpp
	$(CC) $(CXXFLAGS) -c Source/Map.cpp -o Map.o

Atmosphere.o: Source/Atmosphere.cpp
	$(CC) $(CXXFLAGS) -c Source/Atmosphere.cpp -o Atmosphere.o

Earth.o: Source/Earth.cpp Source/Map.cpp Source/Atmosphere.cpp
	$(CC) $(CXXFLAGS) -c Source/Earth.cpp -o Earth.o

SkyCubeMP.o: Source/SkyCubeMP.cpp Source/Earth.cpp Source/Atmosphere.cpp Source/Map.cpp
	$(CC) $(CXXFLAGS) -c Source/SkyCubeMP.cpp -o SkyCubeMP.o

SkyCube.o: Source/SkyCube.cpp Source/Earth.cpp Source/Atmosphere.cpp Source/Map.cpp
	$(CC) $(CXXFLAGS) -c Source/SkyCube.cpp -o SkyCube.o

%.d: %.c
	$(CC) -MM $< > $@

clean: 
	rm -f *.o *~ 
