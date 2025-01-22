#include "Map.h"
#include "Earth.h"
#include "Atmosphere.h"
#include "SkyCube.h"
#include "SkyCubeMP.h"
#include <chrono>

int main()
{
    int nlat = 41;
    // int nlat = 3;
    Earth e = Earth(nlat);
    std::vector<int>& t = e.terrain;
    // std::vector<int> t = std::vector<int>(6*nlat*nlat, 0);
    // t[4] = 1;
    // SkyCube s = SkyCube(nlat, 5, 6378000, 100000, 50000, t);

    // s.Dynamics();

    SkyCubeMP ss = SkyCubeMP(nlat, 4, 6378.000, 100000, 10000, t);
    
    ss.Dynamics(0.1, 1, -1);
    ss.Dynamics(0.1, 6*60*30, 1);

    // auto stop2 = std::chrono::high_resolution_clock::now();
    // auto duration2 = std::chrono::duration_cast<std::chrono::microseconds>(stop2 - start2);

    // std::cout << duration1.count() << std::endl;
    // std::cout << duration2.count() << std::endl;
}
