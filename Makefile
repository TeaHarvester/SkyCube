CC = g++
CXXFLAGS = -Wall -g -fopenmp

ECOSIM: Main.o Map.o Atmosphere.o Earth.o SkyCubeMP.o SkyCube.o
	$(CC) $(CXXFLAGS) Main.o Map.o Atmosphere.o Earth.o SkyCubeMP.o SkyCube.o -o $@

Main.o: Main.cpp Map.cpp Map.h 
	$(CC) $(CXXFLAGS) -c Main.cpp -o Main.o

Map.o: Map.cpp Map.h Atmosphere.cpp Atmosphere.h
	$(CC) $(CXXFLAGS) -c Map.cpp -o Map.o

Atmosphere.o: Atmosphere.cpp Atmosphere.h
	$(CC) $(CXXFLAGS) -c Atmosphere.cpp -o Atmosphere.o

Earth.o: Earth.cpp Earth.h Map.cpp Map.h Atmosphere.cpp Atmosphere.h
	$(CC) $(CXXFLAGS) -c Earth.cpp -o Earth.o

SkyCubeMP.o: SkyCubeMP.cpp SkyCubeMP.h Earth.cpp Earth.h Atmosphere.cpp Atmosphere.h Map.cpp Map.h
	$(CC) $(CXXFLAGS) -c SkyCubeMP.cpp -o SkyCubeMP.o

SkyCube.o: SkyCube.cpp SkyCube.h Earth.cpp Earth.h Atmosphere.cpp Atmosphere.h Map.cpp Map.h
	$(CC) $(CXXFLAGS) -c SkyCube.cpp -o SkyCube.o

clean: 
	rm -f *.o *~ 
