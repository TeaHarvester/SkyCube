#ifndef ATMOSPHERE
#define ATMOSPHERE

// #include <Map>
#include <vector>
#include <random> 
#include <fstream>
#include <string>

class Atmosphere
{
    public:
    void Dynamics();
    Atmosphere(int nlat, int n_lev, double R, double p0, double pt, std::vector<int>& terre);
    Atmosphere(int sx, int sy, int n_lev, double dy_, double p0, double pt, double lt, std::vector<int>& terre);
    ~Atmosphere();

    private:
    double R = 287.05;
    double cp = 1003.5;
    double g = 9.80665;
    double boltz = 5.670/100000000;
    const int xsize;
    const int ysize;
    const int nlevel;
    double dt;
    double ptop;
    double p_ref;
    double OM;
    double radius;
    double ds;
    std::vector<double> dx;
    std::vector<double> dy;
    std::vector<double> ps;
    std::vector<double> ps_buffer;
    std::vector<double> ps_dt;
    std::vector<double> surface_pot;
    std::vector<double> albedo;
    std::vector<double> latitude;
    std::vector<std::vector<double>> T;
    std::vector<std::vector<double>> theta;
    std::vector<std::vector<double>> theta_buffer;
    std::vector<std::vector<double>> theta_dt;
    std::vector<std::vector<double>> geopot;
    std::vector<std::vector<double>> u;
    std::vector<std::vector<double>> v;
    std::vector<std::vector<double>> u_buffer;
    std::vector<std::vector<double>> v_buffer;
    std::vector<std::vector<double>> sigma_vel;
    std::vector<std::vector<double>> u_accel;
    std::vector<std::vector<double>> v_accel;
    std::vector<std::vector<double>> friction;
    std::vector<std::vector<double>> radiant_heatflux;
    std::vector<double> sigma;
    std::vector<int> terrain;
    std::vector<int> inorth;
    std::vector<int> isouth;
    std::vector<int> ieast;
    std::vector<int> iwest;
    std::vector<int> inortheast;
    std::vector<int> inorthwest;
    std::vector<int> isoutheast;
    std::vector<int> isouthwest;
    std::vector<int> ipolar;
    void Hydrostatic();
    void SolvePrimitive();
    void SurfacePressureTendency();
    void SigmaVelocity();
    void ThetaTendency();
    void LeapFrog();
    void EulerForward();
    void Radiate();
    void PrintV(const int level);
    void Print(const int level, char var);
    void PrintP(char var);
    void WriteState();
    std::vector<double> Upstream(int& level, int& pos);
    std::vector<int> GetPeriodicDir(const int centre);
    std::vector<int> GetMercatorDir(const int centre);
};

#endif
