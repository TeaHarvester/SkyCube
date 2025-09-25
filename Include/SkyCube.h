#ifndef SKYCUBE
#define SKYCUBE
#include <fstream>
#include <string>
#include <vector>

class SkyCube
{
    public:
    virtual void Dynamics();
    SkyCube(int nl, int n_lev, double rad, double p0, double pt, std::vector<int>& terre);
    ~SkyCube();

    private:
    double time = 0;
    double R = 287.05;
    double cp = 1003.5;
    double g = 9.80665;
    double OM = 0.000072921;
    double tilt = 7*3.141/180;
    double mass_change = 0;
    int nlat;
    int nlevel;
    double dt;
    double radius;
    double p_ref;
    double ptop;
    double ds;
    double* latitude;
    double* longitude;
    double* sigma;
    double* area;
    double* ps;
    double* ps_buffer;
    double* ps_dt;
    double* surface_pot;
    double* eastlength;
    double* westlength;
    double* northlength;
    double* southlength;
    double* x1i;
    double* x1j;
    double* x2i;
    double* x2j;
    double* y1i;
    double* y1j;
    double* y2i;
    double* y2j;
    double** Wx2i;
    double** Wx2j;  
    double** Ny1i;
    double** Ny1j;
    double* u;
    double* u_buffer;
    double* u_accel;
    double* v;
    double* v_buffer;
    double* v_accel;
    double* geopot;
    double* sigma_vel;
    double* theta;
    double* theta_buffer;
    double* theta_dt;
    double* T;
    double* friction;
    int* inorth;
    int* isouth;
    int* ieast;
    int* iwest;
    int* inortheast;
    int* inorthwest;
    int* isoutheast;
    int* isouthwest;
    double** unorth;
    double** usouth;
    double** ueast;
    double** uwest;
    double** vnorth;
    double** vsouth;
    double** veast;
    double** vwest;
    double** unorthwest;
    double** usouthwest;
    double** vnortheast;
    double** vnorthwest;
    std::vector<int> terrain;
    void SolvePrimitive();
    void Hydrostatic();
    void SurfacePressureTendency();
    void ThetaTendency();
    void Radiate();
    void WriteState();
    void LeapFrog();
    void EulerForward();
    std::vector<double> Upstream(int& level, int& i);
    void Print(int face);
};

#endif