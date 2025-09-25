#ifndef EARTH
#define EARTH
#include"Map.h"

class Earth : public Map
{
    public:
    int nlat;
    Earth(const int nlat_);

    private:
    void Faultline(int tile);
    void Print(int tile);
};


#endif