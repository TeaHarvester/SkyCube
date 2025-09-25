#ifndef MAP
#define MAP

#include <iostream>
#include <vector>
#include <random>
#include <algorithm>
#include <fstream>
#include <cmath>
#include "Atmosphere.h"

class Map
{
    public:
    std::vector<int> terrain;
    Atmosphere* atmos;
    void Print(bool topology);
    Map(const int sx, const int sy);
    ~Map();

    protected:
    const int xsize;
    const int ysize;
    std::vector<int> runoff;
    std::vector<int> snowpack;
    std::vector<int> sources;
    std::vector<int> faultlines; 
    std::vector<char> features;
    std::vector<int> stagger; 
    std::vector<int> Cut(const int centre, bool river, std::vector<int>& upstream);
    std::vector<int> GetDir(const int centre);
    void Mountain(const int centre, const float height, const float width);
    void Glaciate(const int height, bool river);
    void Faultline();
    void DrawTerrain();
    void Write();
    Map(const int nlat);
};

#endif