#include "SkyCube.h"
#include <iostream>
#include <iomanip>
#include <cmath>

void SkyCube::Dynamics()
{
    // WriteState();
    SolvePrimitive();
    EulerForward();
    WriteState();

    double iter = 0;

    // for (int i = 0; i < 24*365; ++i)
    // {
    //     SolvePrimitive();
    //     LeapFrog();
    //     time += dt;
    // }

    while (true)
    {

        // for (int i = 0; i < 60; ++i)
        // {
            SolvePrimitive();
            LeapFrog();
            time += dt;
        // }

        // WriteState();
        ++iter;
    }
}

void SkyCube::SolvePrimitive()
{
    Hydrostatic();
    SurfacePressureTendency();
    ThetaTendency();
    Radiate();

    double k = R/cp;
    int tsize = nlat*nlat;
    int lsize = 6*tsize;

    for (int level = 0; level < nlevel - 1; ++level)
    {
        int l_ = level*lsize;

        for (int i = 0; i < lsize; ++i)
        {
            int pos = i%tsize;
            int epos = ieast[i]%tsize;
            int spos = isouth[i]%tsize;
            double _pi[4];
            double _phi[4];
            double _theta[4];

            _pi[0] = ps[i] - ptop;
            _pi[1] = ps[ieast[i]] - ptop;
            _pi[2] = ps[i] - ptop;
            _pi[3] = ps[isouth[i]] - ptop;
            _phi[0] = geopot[i + l_];
            _phi[1] = geopot[ieast[i] + l_];
            _phi[2] = geopot[i + l_];
            _phi[3] = geopot[isouth[i] + l_];
            _theta[0] = theta[i + l_]/_pi[0];
            _theta[1] = theta[ieast[i] + l_]/_pi[1];
            _theta[2] = theta[i + l_]/_pi[2];
            _theta[3] = theta[isouth[i] + l_]/_pi[3];

            double pix = (_pi[0] + _pi[1])/2;
            double piy = (_pi[2] + _pi[3])/2;

            std::vector<double> advect = Upstream(level, i);

            double sigma_advect[2] = {0, 0};
            double upper_u = level == 0 ? u[i + l_] : u[i + l_ - lsize];
            double upper_v = level == 0 ? v[i + l_] : v[i + l_ - lsize];
            double lower_u = level == nlevel - 2 ? u[i + l_] : u[i + l_ + lsize];
            double lower_v = level == nlevel - 2 ? v[i + l_] : v[i + l_ + lsize];

            double sigma_velx = (sigma_vel[i + l_] + sigma_vel[ieast[i] + l_] + 
                                sigma_vel[i + l_ + lsize] + sigma_vel[ieast[i] + l_ + lsize])/4;
            double sigma_vely = (sigma_vel[i + l_] + sigma_vel[isouth[i] + l_] + 
                                sigma_vel[i + l_ + lsize] + sigma_vel[isouth[i] + l_ + lsize])/4;

            double s_velx_f = (sigma_velx + fabs(sigma_velx))/2;
            double s_velx_b = (sigma_velx - fabs(sigma_velx))/2;
            double s_vely_f = (sigma_vely + fabs(sigma_vely))/2;
            double s_vely_b = (sigma_vely - fabs(sigma_vely))/2;

            sigma_advect[0] = (s_velx_f*(u[i + l_] - lower_u) + s_velx_b*(upper_u - u[i + l_]))/ds;
            sigma_advect[1] = (s_vely_f*(v[i + l_] - lower_v) + s_vely_b*(upper_v - v[i + l_]))/ds;

            double pressure_grad[2];
            double px = (sigma[level] + ds/2)*pix + ptop;
            double p1x = sigma[level]*pix + ptop;
            double p2x = sigma[level + 1]*pix + ptop;
            double py = (sigma[level] + ds/2)*piy + ptop;
            double p1y = sigma[level]*piy + ptop;
            double p2y = sigma[level + 1]*piy + ptop;
            double Tx = (((_theta[0] + _theta[1])/2)/((k + 1)*pow(p_ref, k)*ds*pix))*(pow(p2x, k + 1) - pow(p1x, k + 1));
            double Ty = (((_theta[2] + _theta[3])/2)/((k + 1)*pow(p_ref, k)*ds*piy))*(pow(p2y, k + 1) - pow(p1y, k + 1));
            double rhox = px/(R*Tx);
            double rhoy = py/(R*Ty);

            double& u1i = x1i[pos];
            double& u1j = x1j[pos];
            double& u2i = x2i[pos];
            double& u2j = x2j[pos];
            double& v1i = y1i[pos];
            double& v1j = y1j[pos];
            double& v2i = y2i[pos];
            double& v2j = y2j[pos];
            // double& p1i = pcurvibasis[pos];
            // double& p1j = pcurvibasis[pos + tsize];
            // double& p2i = pcurvibasis[pos + 2*tsize];
            // double& p2j = pcurvibasis[pos + 3*tsize];

            double usin = u1i*u2j - u1j*u2i;
            double vsin = v1i*v2j - v1j*v2i;

            double dpi = u1i*(_pi[1] - _pi[0]);
            //  + u1j*((ps[inorth[i]] + ps[inortheast[i]]) - 
                                                    //   (ps[isouth[i]] + ps[isoutheast[i]]))/2;

            double dpj = v2j*(_pi[2] - _pi[3]);
            //  + v2i*((ps[ieast[i]] + ps[isoutheast[i]]) - 
                                                    //   (ps[iwest[i]] + ps[isouthwest[i]]))/2;

            pressure_grad[0] = (pix/rhox)*(sigma[level] + ds/2)*dpi*usin*(2*eastlength[pos])/(area[pos] + area[epos]);
            pressure_grad[1] = (piy/rhoy)*(sigma[level] + ds/2)*dpj*vsin*(2*southlength[pos])/(area[pos] + area[spos]);

            double geograd[2];
            double dphii = u1i*(_phi[1] - _phi[0]);
            //  + u1j*((geopot[inorth[i] + l_] + geopot[inortheast[i] + l_]) - 
                                                        //   (geopot[isouth[i] + l_] + geopot[isoutheast[i] + l_]))/2;

            double dphij = v2j*(_phi[2] - _phi[3]);
            //  + v2i*((geopot[ieast[i] + l_] + geopot[isoutheast[i] + l_]) - 
                                                        //   (geopot[iwest[i] + l_] + geopot[isouthwest[i] + l_]))/2;

            geograd[0] = pix*dphii*usin*(2*eastlength[pos])/(area[pos] + area[epos]);
            geograd[1] = piy*dphij*vsin*(2*southlength[pos])/(area[pos] + area[spos]);

            u_accel[i + l_] = -advect[0] - sigma_advect[0] - geograd[0] - pressure_grad[0];
            v_accel[i + l_] = -advect[1] - sigma_advect[1] - geograd[1] - pressure_grad[1];
        }
    }
}

void SkyCube::Hydrostatic()
{
    double k = R/cp;
    int lsize = 6*nlat*nlat;

    for (int i = 0; i < lsize; ++i)
    { 
        double pi = ps[i] - ptop;
        double potential = 0;

        for (int level = 0; level < nlevel - 2; ++level)
        {
            int l_ = level*lsize;

            double p0 = pi*sigma[level] + ptop;
            double p = pi*sigma[level + 1] + ptop;
            double p1 = pi*sigma[level + 2] + ptop;
            double PI1 = (cp/((k + 1)*pow(p_ref, k)*pi*ds))*(pow(p, k + 1) - pow(p0, k + 1));
            double PI2 = (cp/((k + 1)*pow(p_ref, k)*pi*ds))*(pow(p1, k + 1) - pow(p, k + 1));
            double PI_dpi = ((cp/(pow(p_ref, k)*ds))*(sigma[level + 1]*pow(p, k) - sigma[level]*pow(p0, k)) - PI1)/pi;
            double theta_inter = (theta[i + l_] + theta[i + l_ + lsize])/(2*pi);

            potential += theta[i + l_]*PI_dpi*ds;
            potential -= sigma[level + 1]*(PI2 - PI1)*theta_inter;
        }

        int level = nlevel - 2;
        int l_ = level*lsize;

        double p0 = pi*sigma[level] + ptop;
        double p = pi*sigma[level + 1] + ptop;
        double PI1 = (cp/((k + 1)*pow(p_ref, k)*pi*ds))*(pow(p, k + 1) - pow(p0, k + 1));
        double PI_dpi = ((cp/(pow(p_ref, k)*ds))*(sigma[level + 1]*pow(p, k) - sigma[level]*pow(p0, k)) - PI1)/pi;

        potential += theta[i + l_]*PI_dpi*ds;
        geopot[i + l_] = surface_pot[i] + potential;

        for (int level = nlevel - 3; level >= 0; --level)
        {
            l_ = level*lsize;
            double p0 = pi*sigma[level] + ptop;
            double p = pi*sigma[level + 1] + ptop;
            double p1 = pi*sigma[level + 2] + ptop;
            double PI1 = (cp/((k + 1)*pow(p_ref, k)*pi*ds))*(pow(p, k + 1) - pow(p0, k + 1));
            double PI2 = (cp/((k + 1)*pow(p_ref, k)*pi*ds))*(pow(p1, k + 1) - pow(p, k + 1));
            double theta_inter = (theta[i + l_] + theta[i + l_ + lsize])/(2*pi);
            geopot[i + l_] = geopot[i + l_ + lsize] + theta_inter*(PI2 - PI1);
        }
    }  
}

void SkyCube::SurfacePressureTendency()
{
    double _u[2];
    double _v[2];
    int tsize = nlat*nlat;
    int lsize = 6*tsize;

    for (int i = 0; i < lsize; ++i)
    {
        int pos = i%tsize;
        // int wpos = iwest[i]%tsize;
        // int npos = inorth[i]%tsize;

        // int nw_corner = !(pos == 0);
        // int sw_corner = !(pos == nlat - 1);
        // int ne_corner = !(pos == nlat*(nlat - 1));
        // int se_corner = !(pos == nlat*nlat - 1);
        int north_halo = 1 - 2*(pos%nlat == 0);
        // int south_halo = 1 - 2*(pos%nlat == nlat - 1);
        int west_halo = 1 - 2*(pos/nlat == 0);

        double pi = ps[i] - ptop;
        // double _u2i = *Wx2i[pos];
        // double& _u2j = *Wx2j[pos];
        // double u2i_ = x2i[pos];
        // double& u2j_ = x2j[pos];
        // double& Nv1i = *Ny1i[pos];
        // double Nv1j = *Ny1j[pos];
        // double& Sv1i = y1i[pos];
        // double Sv1j = y1j[pos];

        double ps_vel[nlevel - 1];
        double cumul_ps_vel[nlevel - 1];
        ps_vel[nlevel - 2] = 0;
        cumul_ps_vel[nlevel - 2] = 0;
        ps_dt[i] = 0;

        for(int level = 0; level < nlevel - 1; ++level)
        {
            int l_ = level*lsize;
            double uw = west_halo*(*(uwest[i] + l_));
            double ue = u[i + l_];
            double vn = north_halo*(*(vnorth[i] + l_));
            double vs = v[i + l_];
            // double vw = (v[i + l_] + west_halo*(*(vwest[i] + l_)) + 
            //             north_halo*((*(vnorth[i] + l_)) + west_halo*(*(vnorthwest[i] + l_))))/4;
            // double ve = (v[i + l_]  + *(veast[i] + l_) +
            //             north_halo*(*(vnorth[i] + l_) + *(vnortheast[i] + l_)))/4;
            // double un = (u[i + l_] + west_halo*(*(uwest[i] + l_)) + 
            //              north_halo*(*(unorth[i] + l_) + west_halo*(*(unorthwest[i] + l_))))/4;
            // double us = (u[i + l_] + west_halo*(*(uwest[i] + l_)) + 
            //              south_halo*(*(usouth[i] + l_) + west_halo*(*(usouthwest[i] + l_))))/4;

            // _u[0] = (_u2j*uw - nw_corner*sw_corner*_u2i*vw);
            // _u[1] = (u2j_*ue - ne_corner*se_corner*u2i_*ve);
            // _v[0] = (Nv1i*vn - nw_corner*ne_corner*Nv1j*un);
            // _v[1] = (Sv1i*vs - sw_corner*se_corner*Sv1j*us);

            _u[0] = uw;
            _u[1] = ue;
            _v[0] = vn;
            _v[1] = vs;

            ps_vel[level] = -ds*(_u[1]*eastlength[pos] - _u[0]*westlength[pos])/area[pos];
            ps_vel[level] -= ds*(_v[0]*northlength[pos] - _v[1]*southlength[pos])/area[pos];
            ps_dt[i] += ps_vel[level];
            cumul_ps_vel[level] = ps_dt[i]; 
            mass_change += ps_vel[level]*area[pos];
        }

        for(int level = 0; level < nlevel - 1; ++level)
        {
            int l_ = level*lsize;
            sigma_vel[i + l_ + lsize] = (sigma[level + 1]*ps_dt[i] - cumul_ps_vel[level])/pi;
        }
    }
}

void SkyCube::ThetaTendency()
{
    // double _theta[4];
    // double _u[2];
    // double _v[2];
    int tsize = nlat*nlat;
    int lsize = 6*nlat*nlat;

    for (int i = 0; i < lsize; ++i)
    {
        int pos = i%tsize;
        // int wpos = iwest[i]%tsize;
        // int npos = inorth[i]%tsize;

        int north_halo = 1 - 2*(pos%nlat == 0);
        // int south_halo = 1 - 2*(pos%nlat == nlat - 1);
        int west_halo = 1 - 2*(pos/nlat == 0);

        double pi = ps[i] - ptop;

        for(int level = 0; level < nlevel - 1; ++level)
        {
            int l_ = level*lsize;
            double uw = west_halo*(*(uwest[i] + l_));
            double ue = u[i + l_];
            double vn = north_halo*(*(vnorth[i] + l_));
            double vs = v[i + l_];

            double uw_f = (uw + fabs(uw))/2;
            double uw_b = (uw - fabs(uw))/2;
            double ue_f = (ue + fabs(ue))/2;
            double ue_b = (ue - fabs(ue))/2;
            double vn_f = (vn + fabs(vn))/2;
            double vn_b = (vn - fabs(vn))/2;
            double vs_f = (vs + fabs(vs))/2;
            double vs_b = (vs - fabs(vs))/2;

            double _theta[4];
            _theta[0] = theta[iwest[i] + l_]/(ps[iwest[i]] - ptop);
            _theta[1] = theta[ieast[i] + l_]/(ps[ieast[i]] - ptop);
            _theta[2] = theta[inorth[i] + l_]/(ps[inorth[i]] - ptop);
            _theta[3] = theta[isouth[i] + l_]/(ps[isouth[i]] - ptop);

            theta_dt[i + l_] = (westlength[pos]*(uw_f*_theta[0] + uw_b*theta[i + l_]/pi) +
                                southlength[pos]*(vs_f*_theta[3] + vs_b*theta[i + l_]/pi) -
                                northlength[pos]*(vn_b*_theta[2] + vn_f*theta[i + l_]/pi) -  
                                eastlength[pos]*(ue_b*_theta[1] + ue_f*theta[i + l_]/pi))/area[pos];

            if (theta[i + l_] + theta_dt[i + l_]*dt < 0)
            {
                int stop = 1;
            }

        }

        for (int level = 0; level < nlevel - 2; ++level)
        {
            int l_ = level*lsize;
            double s_vel = sigma_vel[i + l_ + lsize];
            double s_vel_f = (s_vel + fabs(s_vel))/2;
            double s_vel_b = (s_vel - fabs(s_vel))/2;
            double& upper_theta = theta[i + l_];
            double& lower_theta = theta[i + l_ + lsize];
            double sigma_advect_upper = (s_vel_f*lower_theta + s_vel_b*upper_theta)/ds;
            double sigma_advect_lower = -sigma_advect_upper;
            theta_dt[i + l_] = theta_dt[i + l_] + sigma_advect_upper;
            theta_dt[i + l_ + lsize] = theta_dt[i + l_ + lsize] + sigma_advect_lower;

            if (theta[i + l_] + theta_dt[i + l_]*dt < 0 || i == 132)
            {
                int stop = 1;
            }
        }
    }
}

void SkyCube::Radiate()
{
    int tsize = nlat*nlat;
    int lsize = 6*tsize;

    for (int i = 0; i < lsize; ++i)
    {
        int pos = i%tsize;
        double pi = ps[i] - ptop;
        double irradiance = 1362;
        double em = 0.78;
        double boltz = 5.670/100000000;
        double azimuth = OM*time;
        double tilt_ = 0.4*sin(OM*time/(365*24*60*60));
        double albedo = 0.1;
        double GH = 0.5;
        double T_to_theta = theta[i]/(T[i]*pi);
        albedo -= (terrain[i] > 0)*0.025;
        albedo += (T[i + (nlevel - 2)*lsize] < 263)*0.2;
        albedo = 2*((terrain[i] > 0) && (T[i + (nlevel - 2)*lsize] < 263))*albedo + 
                    !(((terrain[i] > 0) && (T[i + (nlevel - 2)*lsize] < 263)))*albedo;
        double irrad = cos(longitude[i] + azimuth) > 0; 
        double rad_flux = dt*irrad*(1 - albedo)*irradiance*T_to_theta*(latitude[i] - tilt_)*std::fabs(cos(longitude[i] + azimuth))*g/cp;

        theta_dt[i + (nlevel - 2)*lsize] += rad_flux;

        for (int level = 0; level < nlevel - 1; ++level)
        {
            int l_ = level*lsize;

            if (theta[i + l_] + theta_dt[i + l_]*dt < 0)
            {
                int stop = 1;
            }

            theta_dt[i + l_] -= dt*T_to_theta*GH*em*boltz*pow(T[i + l_], 4)*g/cp;
            double loss = dt*T_to_theta*GH*em*boltz*pow(T[i + l_], 4)*g/cp;

            if (theta[i + l_] + theta_dt[i + l_]*dt < 0)
            {
                int stop = 1;
            }
        }
    }
}

std::vector<double> SkyCube::Upstream(int& level, int& i)
{
    int tsize = nlat*nlat;
    int pos = i%tsize;
    // int wpos = iwest[i]%tsize;
    // int npos = inorth[i]%tsize;
    int epos = ieast[i]%tsize;
    int spos = isouth[i]%tsize;

    int north_halo = 1 - 2*(pos%nlat == 0);
    int south_halo = 1 - 2*(pos%nlat == nlat - 1);
    int east_halo = 1 - 2*(pos/nlat == nlat - 1);
    int west_halo = 1 - 2*(pos/nlat == 0);
    int northeast_halo = (pos != nlat*(nlat - 1))*north_halo;
    int southwest_halo = (pos != nlat - 1)*south_halo;

    double _pi[4];

    // double& u1i = x1i[i];
    // double& u1j = x1j[i];
    // double& u2i = x2i[i];
    // double& u2j = x2j[i];
    // double& v1i = y1i[i];
    // double& v1j = y1j[i];
    // double& v2i = y2i[i];
    // double& v2j = y2j[i];

    _pi[0] = ps[i] - ptop;
    _pi[1] = ps[ieast[i]] - ptop;
    _pi[2] = ps[i] - ptop;
    _pi[3] = ps[isouth[i]] - ptop;

    int lsize = 6*tsize;
    int l_ = level*lsize;

    double uw = west_halo*(*(uwest[i] + l_));
    double ue = east_halo*(*(ueast[i] + l_));
    double un = *(unorth[i] + l_);
    double us = south_halo*(*(usouth[i] + l_));
    double vw = west_halo*(*(vwest[i] + l_));
    double ve = *(veast[i] + l_);
    double vn = north_halo*(*(vnorth[i] + l_));
    double vs = south_halo*(*(vsouth[i] + l_));

    double upper_v = (vn + northeast_halo*(*(vnortheast[i] + l_)))/2;
    double lower_v = (v[i + l_] + ve)/2;
    double west_u = (uw + southwest_halo*(*(usouthwest[i] + l_)))/2;
    double east_u = (u[i + l_] + us)/2;

    double upper_pi = (ps[ieast[i]] + ps[i] + ps[inortheast[i]] + ps[inorth[i]] - 4*ptop)/4;
    double lower_pi = (ps[ieast[i]] + ps[i] + ps[isoutheast[i]] + ps[isouth[i]] - 4*ptop)/4;
    double west_pi = (ps[iwest[i]] + ps[i] + ps[isouthwest[i]] + ps[isouth[i]] - 4*ptop)/4;
    double east_pi = (ps[ieast[i]] + ps[i] + ps[isoutheast[i]] + ps[isouth[i]] - 4*ptop)/4;

    // double ux = u2j*(2*u[i + l_] + ue + uw)/4 - u2i*(upper_v + lower_v);
    // double vy = v1i*(2*v[i + l_] + vn + vs)/4 - v1j*(east_u + west_u);

    // double ux = (2*u[i + l_] + ue + uw)/4;
    // double vy = (2*v[i + l_] + vn + vs)/4;
    double ux = (2*u[i + l_] + ue + uw)/2;
    double vy = (2*v[i + l_] + vn + vs)/2;

    // upper_v = (u1i*upper_v - u1j*(u[i + l_] + un)/2);
    // lower_v = (u1i*lower_v - u1j*(u[i + l_] + us)/2);
    // west_u = (v2j*west_u - v2i*(v[i + l_] + vw)/2);
    // east_u = (v2j*east_u - v2i*(v[i + l_] + ve)/2);

    double ux_f = (ux + fabs(ux))/(_pi[0] + _pi[1]);
    double ux_b = (ux - fabs(ux))/(_pi[0] + _pi[1]);
    double vx_f = (upper_v + lower_v + fabs(upper_v + lower_v))/(2*(upper_pi + lower_pi));
    double vx_b = (upper_v + lower_v - fabs(upper_v + lower_v))/(2*(upper_pi + lower_pi));

    double vy_f = (vy + fabs(vy))/(_pi[2] + _pi[3]);
    double vy_b = (vy - fabs(vy))/(_pi[2] + _pi[3]);
    double uy_f = (west_u + east_u + fabs(west_u + east_u))/(2*(west_pi + east_pi));
    double uy_b = (west_u + east_u - fabs(west_u + east_u))/(2*(west_pi + east_pi));

    double u_advect = 2*(eastlength[pos]*(ux_f*(u[i + l_] - uw) + ux_b*(ue - u[i + l_])) + 
                         (southlength[pos]*vx_f*(u[i + l_] - us) + northlength[pos]*vx_b*(un - u[i + l_])))
                         /(area[pos] + area[epos]);

    double v_advect = 2*(southlength[pos]*(vy_f*(v[i + l_] - vs) + vy_b*(vn - v[i + l_])) + 
                         (westlength[pos]*uy_f*(v[i + l_] - vw) + eastlength[pos]*uy_b*(ve - v[i + l_])))
                         /(area[pos] + area[spos]);

    std::vector<double> advect = {u_advect, v_advect};

    return advect;
}

void SkyCube::LeapFrog()
{
    int tsize = nlat*nlat;
    int lsize = 6*tsize;
    double p_swap;
    double u_swap;
    double v_swap;
    double theta_swap;
    double timefilter = 0.01;

    for (int i = 0; i < lsize; ++i)
    {
        p_swap = ps[i];
        ps[i] = ps_buffer[i] + ps_dt[i]*dt; 
        ps_buffer[i] = p_swap + timefilter*(ps[i] -2*p_swap + ps_buffer[i]);  

        if (ps[i] < 80000)
        {
            bool stop = true;
        }
    }

    for (int level = 0; level < nlevel - 1; ++level)
    {
        int l_ = level*lsize;

        for (int i = 0; i < lsize; ++i)
        {
            int pos = i%tsize;
            int north_halo = 1 - 2*(pos%nlat == 0);
            int south_halo = 1 - 2*(pos%nlat == nlat - 1);
            int east_halo = 1 - 2*(pos/nlat == nlat - 1);
            int west_halo = 1 - 2*(pos/nlat == 0);
            int northeast_halo = (pos != nlat*(nlat - 1))*north_halo;
            int southwest_halo = (pos != nlat - 1)*south_halo;

            u_swap = u[i + l_];
            v_swap = v[i + l_];
            theta_swap = theta[i + l_];
            u[i + l_] = u_buffer[i + l_] + u_accel[i + l_]*dt;
            v[i + l_] = v_buffer[i + l_] + v_accel[i + l_]*dt;

            if (theta[i + l_] + theta_dt[i + l_]*dt < 0)
            {
                int stop = 1;
            }

            double pi = ps[i] - ptop;
            theta[i + l_] = theta_buffer[i + l_] + theta_dt[i + l_]*dt;
            T[i + l_] = theta[i + l_]*pow((p_ref/((sigma[level] + ds/2)*(ps[i] - ptop) + ptop)), -R/cp)/pi;

            // implicit schemes
            double f = 2*OM*sin(latitude[i]);

            double u_interp = (u[i + l_] + west_halo*(*(uwest[i] + l_))
             + south_halo*(*(usouth[i] + l_)) + southwest_halo*(*(usouthwest[i] + l_)))/4;

            double v_interp = (v[i + l_] + north_halo*(*(vnorth[i] + l_))
             + east_halo*(*(veast[i] + l_)) + northeast_halo*(*(vnortheast[i] + l_)))/4;

            u[i + l_] = (u[i + l_] - f*dt*v_interp)/(pow(f*dt, 2) + friction[i + l_] + 1);
            v[i + l_] = (v[i + l_] + f*dt*u_interp)/(pow(f*dt, 2) + friction[i + l_] + 1);

            u_buffer[i + l_] = u_swap + timefilter*(u[i + l_] -2*u_swap + u_buffer[i + l_]);
            v_buffer[i + l_] = v_swap + timefilter*(v[i + l_] -2*v_swap + v_buffer[i + l_]);
            theta_buffer[i + l_] = theta_swap + timefilter*(theta[i + l_] -2*theta_swap + theta_buffer[i + l_]);
        }
    }
}

void SkyCube::EulerForward()
{
    int lsize = 6*nlat*nlat;

    for (int i = 0; i < lsize; ++i)
    {
        ps[i] = ps[i] + ps_dt[i]*dt/2; 
    }

    for (int level = 0; level < nlevel - 1; ++level)
    {
        int l_ = level*lsize;

        for (int i = 0; i < lsize; ++i)
        {
            u[i + l_] = u[i + l_] + u_accel[i + l_]*dt/2;
            v[i + l_] = v[i + l_] + v_accel[i + l_]*dt/2;

            double pi = ps[i] - ptop;
            theta[i + l_] = theta[i + l_] + theta_dt[i + l_]*dt/2;
            T[i + l_] = theta[i + l_]*pow((p_ref/((sigma[level] + ds/2)*(ps[i] - ptop) + ptop)), -R/cp)/pi;

            if (theta[i + l_] < 0)
            {
                int stop = 1;
            }
        }   
    }

    // SurfacePressureTendency();
    // ThetaTendency();
    

    // for (int level = 0; level < nlevel - 1; ++level)
    // {
    //     int l_ = level*lsize;

    //     for (int i = 0; i < lsize; ++i)
    //     {
    //         double pi = ps[i] - ptop;
    //         theta[i + l_] = theta[i + l_] + theta_dt[i + l_]*dt/2;
    //         T[i + l_] = theta[i + l_]*pow((p_ref/((sigma[level] + ds/2)*(ps[i] - ptop) + ptop)), -R/cp)/pi;
    //     }
    // }
}

void SkyCube::Print(int face)
{
    int tsize = nlat*nlat;

    std::cout << std::fixed;
    std::cout << std::endl;

    for (int i = 0; i < 2*nlat; ++i)
    {
        std::cout << std::endl;

        for (int j = 0; j < nlat; ++j)
        {
            int pos = i*nlat/2 + j;

            if (i%2 == 0)
            {
                std::cout << std::scientific << std::setprecision(12) << ps_dt[i/2 + j*nlat + face*tsize]*area[pos] << " ";
                //  << " " << std::setprecision(2) << u[i/2 + j*nlat + face*tsize]*eastlength[pos] << " ";
            }

            // else
            // {
            //     std::cout << std::setprecision(2) << v[i/2 + j*nlat + face*tsize]*southlength[pos] << "      ";
            // }
        }
    }
}

void SkyCube::WriteState()
{
    std::ofstream ufile("u.csv", std::ios_base::app);
    std::ofstream vfile("v.csv", std::ios_base::app);
    std::ofstream sfile("s.csv", std::ios_base::app);
    std::ofstream pfile("p.csv", std::ios_base::app);
    std::ofstream gfile("g.csv", std::ios_base::app);
    std::ofstream tfile("t.csv", std::ios_base::app);
    std::ostringstream streamObj;
    std::string strObj;
    
    int tsize = nlat*nlat;
    int lsize = 6*tsize;

    for (int level = 0; level < nlevel - 1; ++level)
    {
        int l_ = level*lsize;

        for (int i = 0; i < lsize; ++i)
        {
            int pos = i%tsize;
            int pi = ps[i] - ptop;
            int north_halo = 1 - 2*(pos%nlat == 0);
            int west_halo = 1 - 2*(pos/nlat == 0);

            streamObj << std::scientific << std::setprecision(4);
            streamObj << (u[i + l_] + west_halo*(*(uwest[i] + l_)))/(2*(ps[i] - ptop));
            strObj = streamObj.str();
            ufile << strObj << ',';
            streamObj.str(std::string());

            streamObj << (v[i + l_] + north_halo*(*(vnorth[i] + l_)))/(2*(ps[i] - ptop));
            strObj = streamObj.str();
            vfile << strObj << ',';
            streamObj.str(std::string());

            streamObj << std::scientific << std::setprecision(4);
            streamObj << sigma_vel[i + l_];
            strObj = streamObj.str();
            sfile << strObj << ',';
            streamObj.str(std::string());

            streamObj << std::scientific << std::setprecision(7);
            streamObj << geopot[i + l_];
            strObj = streamObj.str();
            gfile << strObj << ',';
            streamObj.str(std::string());

            streamObj << std::scientific << std::setprecision(5);
            streamObj << theta[i + l_]/pi;
            strObj = streamObj.str();
            tfile << strObj << ',';
            streamObj.str(std::string());
        }

        ufile << std::endl;
        vfile << std::endl;
        sfile << std::endl;
        gfile << std::endl;
        tfile << std::endl;
    }

    streamObj << std::scientific << std::setprecision(6);

    for (int i = 0; i < lsize; ++i)
    {
        streamObj << ps[i];
        strObj = streamObj.str();
        pfile << strObj << ',';
        streamObj.str(std::string()); 
    }

    ufile << std::endl;
    vfile << std::endl;
    sfile << std::endl;
    gfile << std::endl;
    pfile << std::endl;
    tfile << std::endl;

    ufile.close();
    vfile.close();
    sfile.close();
    gfile.close();
    pfile.close();
    tfile.close();
}

SkyCube::SkyCube(int nl, int n_lev, double rad, double p0, double pt, std::vector<int>& terre)
:
nlat(nl),
nlevel(n_lev),
dt(60),
radius(rad),
p_ref(p0),
ptop(pt),
ds(1/((double)nlevel - 1)),
latitude(new double[6*nlat*nlat]),
longitude(new double[6*nlat*nlat]),
sigma(new double[nlevel]),
area(new double[nlat*nlat]),
ps(new double[6*nlat*nlat]),
ps_buffer(new double[6*nlat*nlat]),
ps_dt(new double[6*nlat*nlat]),
surface_pot(new double[6*nlat*nlat]),
eastlength(new double[nlat*nlat]),
westlength(new double[nlat*nlat]),
northlength(new double[nlat*nlat]),
southlength(new double[nlat*nlat]),
x1i(new double[6*nlat*nlat]),
x1j(new double[6*nlat*nlat]),
x2i(new double[6*nlat*nlat]),
x2j(new double[6*nlat*nlat]),
y1i(new double[6*nlat*nlat]),
y1j(new double[6*nlat*nlat]),
y2i(new double[6*nlat*nlat]),
y2j(new double[6*nlat*nlat]),
Wx2i(new double*[6*nlat*nlat]),
Wx2j(new double*[6*nlat*nlat]),  
Ny1i(new double*[6*nlat*nlat]),
Ny1j(new double*[6*nlat*nlat]),
u(new double[6*nlat*nlat*(nlevel - 1)]),
u_buffer(new double[6*nlat*nlat*(nlevel - 1)]),
u_accel(new double[6*nlat*nlat*(nlevel - 1)]),
v(new double[6*nlat*nlat*(nlevel - 1)]),
v_buffer(new double[6*nlat*nlat*(nlevel - 1)]),
v_accel(new double[6*nlat*nlat*(nlevel - 1)]),
geopot(new double[6*nlat*nlat*(nlevel - 1)]),
sigma_vel(new double[6*nlat*nlat*nlevel]),
theta(new double[6*nlat*nlat*(nlevel - 1)]),
theta_buffer(new double[6*nlat*nlat*(nlevel - 1)]),
theta_dt(new double[6*nlat*nlat*(nlevel - 1)]),
T(new double[6*nlat*nlat*(nlevel - 1)]),
friction(new double[6*nlat*nlat*(nlevel - 1)]),
inorth(new int[6*nlat*nlat]),
isouth(new int[6*nlat*nlat]),
ieast(new int[6*nlat*nlat]),
iwest(new int[6*nlat*nlat]),
inortheast(new int[6*nlat*nlat]),
inorthwest(new int[6*nlat*nlat]),
isoutheast(new int[6*nlat*nlat]),
isouthwest(new int[6*nlat*nlat]),
unorth(new double*[6*nlat*nlat]),
usouth(new double*[6*nlat*nlat]),
ueast(new double*[6*nlat*nlat]),
uwest(new double*[6*nlat*nlat]),
vnorth(new double*[6*nlat*nlat]),
vsouth(new double*[6*nlat*nlat]),
veast(new double*[6*nlat*nlat]),
vwest(new double*[6*nlat*nlat]),
unorthwest(new double*[6*nlat*nlat]),
usouthwest(new double*[6*nlat*nlat]),
vnortheast(new double*[6*nlat*nlat]),
vnorthwest(new double*[6*nlat*nlat]),
terrain(terre)
{
    int tsize = nlat*nlat;
    int lsize = 6*tsize;
    double pi = 3.1415926535;
    double alpha_incr = pi/(2*nlat);
    double a = 1;

    for (int i = 0; i < nlat; ++i)
    {
        double alpha = pi*(1 - nlat)/(4*nlat) + i*alpha_incr;
        double alpha_edge1 = -pi/4 + i*alpha_incr;
        double alpha_edge2 = -pi/4 + (i + 1)*alpha_incr;
        double x = tan(alpha);
        double _x = tan(alpha_edge1);
        double x_ = tan(alpha_edge2);

        for (int j = 0; j < nlat; ++j)
        {
            int pos = i*nlat + j;

            double beta = pi*(1 - nlat)/(4*nlat) + (nlat - j - 1)*alpha_incr;
            double beta_edge1 = -pi/4 + (nlat - j - 1)*alpha_incr;
            double beta_edge2 = -pi/4 + (nlat - j)*alpha_incr;
            double y = tan(beta);
            double y1 = tan(beta_edge1);
            double y2 = tan(beta_edge2);

            double r = sqrt(x*x + y*y + a*a);
            double _r = sqrt(_x*_x + y*y + a*a);
            double r_ = sqrt(x_*x_ + y*y + a*a);
            double r1 = sqrt(x*x + y1*y1 + a*a);
            double r2 = sqrt(x*x + y2*y2 + a*a);
            double _r1 = sqrt(_x*_x + y1*y1 + a*a);
            double _r2 = sqrt(_x*_x + y2*y2 + a*a);
            double r_1 = sqrt(x_*x_ + y1*y1 + a*a);
            double r_2 = sqrt(x_*x_ + y2*y2 + a*a);

            double e1i = (x_*radius/r_ - _x*radius/_r)/cos(atan((a*radius/r_ - a*radius/_r)/(x_*radius/r_ - _x*radius/_r)));
            double e1j = (y*radius/r_ - y*radius/_r)/cos(atan((a*radius/r_ - a*radius/_r)/(y*radius/r_ - y*radius/_r)));
            double e2i = (x*radius/r2 - x*radius/r1)/cos(atan((a*radius/r2 - a*radius/r1)/(x*radius/r2 - x*radius/r1)));
            double e2j = (y2*radius/r2 - y1*radius/r1)/cos(atan((a*radius/r2 - a*radius/r1)/(y2*radius/r2 - y1*radius/r1)));

            e1i = e1i != e1i ? 0 : e1i;
            e1j = e1j != e1j ? 0 : e1j;
            e2i = e2i != e2i ? 0 : e2i;
            e2j = e2j != e2j ? 0 : e2j;

            area[pos] = e1i*e2j - e1j*e2i;

            // pcurvibasis[pos] = e1i/sqrt(e1i*e1i + e1j*e1j);
            // pcurvibasis[pos + tsize] = e1j/sqrt(e1i*e1i + e1j*e1j);
            // pcurvibasis[pos + 2*tsize] = e2i/sqrt(e2i*e2i + e2j*e2j);
            // pcurvibasis[pos + 3*tsize] = e2j/sqrt(e2i*e2i + e2j*e2j);

            eastlength[pos] = sqrt(pow(x_*radius/r_2 - x_*radius/r_1, 2) + 
                                          pow(y2*radius/r_2 - y1*radius/r_1, 2) + 
                                          pow(a*radius/r_2 - a*radius/r_1, 2));
            westlength[pos] = sqrt(pow(_x*radius/_r2 - _x*radius/_r1, 2) + 
                                          pow(y2*radius/_r2 - y1*radius/_r1, 2) + 
                                          pow(a*radius/_r2 - a*radius/_r1, 2));
            northlength[pos] = sqrt(pow(x_*radius/r_2 - _x*radius/_r2, 2) + 
                                           pow(y2*radius/r_2 - y2*radius/_r2, 2) + 
                                           pow(a*radius/r_2 - a*radius/_r2, 2));
            southlength[pos] = sqrt(pow(x_*radius/r_1 - _x*radius/_r1, 2) + 
                                           pow(y1*radius/r_1 - y1*radius/_r1, 2) + 
                                           pow(a*radius/r_1 - a*radius/_r1, 2));

            double x_e = tan(pi*(1 - nlat)/(4*nlat) + ((i + 1)%nlat)*alpha_incr);
            double r_e = sqrt(x_e*x_e + y*y + a*a);
            double y_s = tan(pi*(1 - nlat)/(4*nlat) + ((nlat - j - 2)%nlat)*alpha_incr);
            double r_s = sqrt(x*x + y_s*y_s + a*a);

            double e1i_ = (x_e*radius/r_e - x*radius/r)/cos(atan((a*radius/r_e - a*radius/r)/(x_e*radius/r_e - x*radius/r)));
            double e1j_ = (y*radius/r_e - y*radius/r)/cos(atan((a*radius/r_e - a*radius/r)/(y*radius/r_e - y*radius/r)));
            double e2i_ = (x*radius/r - x*radius/r_s)/cos(atan((a*radius/r - a*radius/r_s)/(x*radius/r - x*radius/r_s)));
            double e2j_ = (y*radius/r - y_s*radius/r_s)/cos(atan((a*radius/r - a*radius/r_s)/(y*radius/r - y_s*radius/r_s)));
            
            e1i_ = e1i_ != e1i_ ? 0 : e1i_;
            e1j_ = e1j_ != e1j_ ? 0 : e1j_;
            e2i_ = e2i_ != e2i_ ? 0 : e2i_;
            e2j_ = e2j_ != e2j_ ? 0 : e2j_;

            e1i = x_*radius/r_1 - _x*radius/_r1;
            e1j = y1*radius/r_1 - y1*radius/_r1;
            e2i = x_*radius/r_2 - x_*radius/r_1;
            e2j = y2*radius/r_2 - y1*radius/r_1;

            x1i[pos] = e1i_/sqrt(e1i_*e1i_ + e1j_*e1j_);
            x1j[pos] = e1j_/sqrt(e1i_*e1i_ + e1j_*e1j_);
            x2i[pos] = e2i/sqrt(e2i*e2i + e2j*e2j);
            x2j[pos] = e2j/sqrt(e2i*e2i + e2j*e2j);
            y1i[pos] = e1i/sqrt(e1i*e1i + e1j*e1j);
            y1j[pos] = e1j/sqrt(e1i*e1i + e1j*e1j);
            y2i[pos] = e2i_/sqrt(e2i_*e2i_ + e2j_*e2j_);
            y2j[pos] = e2j_/sqrt(e2i_*e2i_ + e2j_*e2j_);

            if (i == nlat - 1)
            {
                x1i[pos] = 1;
                x1j[pos] = 0;
            }

            if (j == nlat - 1)
            {
                y2i[pos] = 0;
                y2j[pos] = 1;
            }

            // x1i[pos + 2*tsize] = x2j[pos];
            // x1j[pos + 2*tsize] = x2i[pos];
            // x2i[pos + 2*tsize] = x1j[pos];
            // x2j[pos + 2*tsize] = x1i[pos];
            // y1i[pos + 2*tsize] = y2j[pos];
            // y1j[pos + 2*tsize] = y2i[pos];
            // y2i[pos + 2*tsize] = y1j[pos];
            // y2j[pos + 2*tsize] = y1i[pos];

            longitude[pos] = alpha;
            longitude[pos + tsize] = pi/2 - beta;
            longitude[pos + 2*tsize] = alpha + 2*pi/2;
            longitude[pos + 3*tsize] =  3*pi/2 - beta;
            longitude[pos + 4*tsize] = -atan(tan(alpha)/tan(beta)) + (j > nlat/2 ? 3*pi/2 : pi/2);
            longitude[pos + 5*tsize] = -atan(tan(-alpha)/tan(-beta)) + (j > nlat/2 ? pi/2 : 3*pi/2);

            latitude[pos] = atan(tan(beta)*cos(alpha));
            latitude[pos + tsize] = atan(tan(alpha)*cos(beta));
            latitude[pos + 2*tsize] = -atan(tan(beta)*cos(alpha));
            latitude[pos + 3*tsize] = -atan(tan(alpha)*cos(beta));
            latitude[pos + 4*tsize] = atan(cos(longitude[pos + 4*tsize])/tan(alpha));
            latitude[pos + 5*tsize] = atan(cos(longitude[pos + 5*tsize])/tan(alpha));

            if (j == nlat/2)
            {
                // longitude[pos + 4*tsize] += pi;
                longitude[pos + 5*tsize] = 2*pi - longitude[pos + 5*tsize];
                latitude[pos + 4*tsize] = atan(cos(longitude[pos + 4*tsize])/tan(alpha));
                latitude[pos + 5*tsize] = atan(cos(longitude[pos + 5*tsize])/tan(alpha));
            }

            if (i == nlat/2)
            {
                latitude[pos + 4*tsize] = latitude[j*nlat + i + 4*tsize];
                latitude[pos + 5*tsize] = latitude[j*nlat + i + 5*tsize];
            }
        }  
    }

    if (nlat%2)
    {
        latitude[nlat*(nlat/2) + nlat/2 + 4*tsize] = pi/2;
        latitude[nlat*(nlat/2) + nlat/2 + 5*tsize] = -pi/2;
        longitude[nlat*(nlat/2) + nlat/2 + 4*tsize] = 0;
        longitude[nlat*(nlat/2) + nlat/2 + 5*tsize] = 0;
    }

    for (int j = 0; j < nlat; ++j)
    {
        latitude[j + nlat*(nlat/2) + 4*tsize] = latitude[j*nlat + nlat/2 + 4*tsize];
        latitude[j + nlat*(nlat/2) + 5*tsize] = latitude[j*nlat + nlat/2 + 5*tsize];
    }

    for (int level = 0; level < nlevel - 1; ++level)
    {
        int l_ = lsize*level;;

        for (int i = 0; i < lsize; ++i)
        {
            u[i + l_] = 0;
            u_buffer[i + l_] = 0;
            u_accel[i + l_] = 0;
            v[i + l_] = 0;
            v_buffer[i + l_] = 0;
            v_accel[i + l_] = 0;
            sigma_vel[i + l_] = 0;
            theta[i + l_] = 0;
            theta_buffer[i + l_] = 0;
            theta_dt[i + l_] = 0;
        }
    }

    double Ts = 273;

    for (int i = 0; i < lsize; ++i)
    {
        int l_ = lsize*(nlevel - 1);

        sigma_vel[i + l_] = 0;
        surface_pot[i] = (R/100000)*terrain[i]*g;
        ps[i] = pow(pow(p_ref, R/cp)*(1 - surface_pot[i]/(cp*Ts)), cp/R);
        ps_buffer[i] = ps[i];
        ps_dt[i] = 0;
    }

    for (int level = 0; level < nlevel; ++level)
    {
        sigma[level] = (double)level/((double)nlevel - 1);
    }

    for (int i = 0; i < nlat; ++i)
    {
        for (int j = 0; j < nlat; ++j)
        {
            int pos = i*nlat + j;
            int tsize = nlat*nlat;

            for (int p = 0; p < 6; ++p)
            {
                if (i > 0)
                {
                    iwest[pos + p*tsize] = pos - nlat + p*tsize;
                    uwest[pos + p*tsize] = &u[iwest[pos + p*tsize]];
                    vwest[pos + p*tsize] = &v[iwest[pos + p*tsize]];

                    Wx2i[pos + p*tsize] = &x2i[iwest[pos + p*tsize]];
                    Wx2j[pos + p*tsize] = &x2j[iwest[pos + p*tsize]];
                }

                if (i < nlat)
                {
                    ieast[pos + p*tsize] = pos + nlat + p*tsize;
                    ueast[pos + p*tsize] = &u[ieast[pos + p*tsize]];
                    veast[pos + p*tsize] = &v[ieast[pos + p*tsize]];
                }

                if (j > 0)
                {
                    inorth[pos + p*tsize] = pos - 1 + p*tsize;
                    unorth[pos + p*tsize] = &u[inorth[pos + p*tsize]];
                    vnorth[pos + p*tsize] = &v[inorth[pos + p*tsize]];

                    Ny1i[pos + p*tsize] = &y1i[inorth[pos + p*tsize]];
                    Ny1j[pos + p*tsize] = &y1j[inorth[pos + p*tsize]];
                }

                if (j < nlat)
                {
                    isouth[pos + p*tsize] = pos + 1 + p*tsize;
                    usouth[pos + p*tsize] = &u[isouth[pos + p*tsize]];
                    vsouth[pos + p*tsize] = &v[isouth[pos + p*tsize]];
                }
            }

            if (i == nlat - 1)
            {
                ieast[pos] = tsize + (nlat - j - 1)*nlat;
                ieast[pos + tsize] = 4*tsize + (nlat - j - 1)*nlat;
                ieast[pos + 2*tsize] = 3*tsize + (nlat - j - 1)*nlat;
                ieast[pos + 3*tsize] = 5*tsize + (nlat - j - 1)*nlat;
                ieast[pos + 4*tsize] = (nlat - j - 1)*nlat;
                ieast[pos + 5*tsize] = 2*tsize + (nlat - j - 1)*nlat;

                ueast[pos] = &v[ieast[pos]];
                ueast[pos + tsize] = &v[ieast[pos + tsize]];
                ueast[pos + 2*tsize] = &v[ieast[pos + 2*tsize]];
                ueast[pos + 3*tsize] = &v[ieast[pos + 3*tsize]];
                ueast[pos + 4*tsize] = &v[ieast[pos + 4*tsize]];
                ueast[pos + 5*tsize] = &v[ieast[pos + 5*tsize]];

                veast[pos] = &u[iwest[ieast[pos]]];
                veast[pos + tsize] = &u[iwest[ieast[pos + tsize]]];
                veast[pos + 2*tsize] = &u[iwest[ieast[pos + 2*tsize]]];
                veast[pos + 3*tsize] = &u[iwest[ieast[pos + 3*tsize]]];
                veast[pos + 4*tsize] = &u[iwest[ieast[pos + 4*tsize]]];
                veast[pos + 5*tsize] = &u[iwest[ieast[pos + 5*tsize]]];

                if (j == nlat - 1)
                {
                    veast[pos] = &v[nlat - 1 + 5*tsize];
                    veast[pos + tsize] = &v[nlat - 1 + 2*tsize];
                    veast[pos + 2*tsize] = &v[nlat - 1 + 4*tsize];
                    veast[pos + 3*tsize] = &v[nlat - 1];
                    veast[pos + 4*tsize] = &v[nlat - 1 + 3*tsize];
                    veast[pos + 5*tsize] = &v[nlat - 1 + tsize];
                }
            }

            if (i == 0)
            {
                iwest[pos] = 3*tsize + (j + 1)*nlat - 1;
                iwest[pos + tsize] = 5*tsize + (j + 1)*nlat - 1;
                iwest[pos + 2*tsize] = tsize + (j + 1)*nlat - 1;
                iwest[pos + 3*tsize] = 4*tsize + (j + 1)*nlat - 1;
                iwest[pos + 4*tsize] = 2*tsize + (j + 1)*nlat - 1;
                iwest[pos + 5*tsize] = (j + 1)*nlat - 1;

                uwest[pos] = &v[iwest[pos]];
                uwest[pos + tsize] = &v[iwest[pos + tsize]];
                uwest[pos + 2*tsize] = &v[iwest[pos + 2*tsize]];
                uwest[pos + 3*tsize] = &v[iwest[pos + 3*tsize]];
                uwest[pos + 4*tsize] = &v[iwest[pos + 4*tsize]];
                uwest[pos + 5*tsize] = &v[iwest[pos + 5*tsize]];

                vwest[pos] = &u[iwest[pos]];
                vwest[pos + tsize] = &u[iwest[pos + tsize]];
                vwest[pos + 2*tsize] = &u[iwest[pos + 2*tsize]];
                vwest[pos + 3*tsize] = &u[iwest[pos + 3*tsize]];
                vwest[pos + 4*tsize] = &u[iwest[pos + 4*tsize]];
                vwest[pos + 5*tsize] = &u[iwest[pos + 5*tsize]];

                Wx2i[pos] = &y1j[iwest[pos]%tsize];
                // Wx2i[pos + tsize] = &y2i[iwest[pos + tsize]];
                // Wx2i[pos + 2*tsize] = &y2i[iwest[pos + 2*tsize]];
                // Wx2i[pos + 3*tsize] = &y2i[iwest[pos + 3*tsize]];
                // Wx2i[pos + 4*tsize] = &y2i[iwest[pos + 4*tsize]];
                // Wx2i[pos + 5*tsize] = &y2i[iwest[pos + 5*tsize]];
                Wx2j[pos] = &y1i[iwest[pos]%tsize];
                // Wx2j[pos + tsize] = &y2j[iwest[pos + tsize]];
                // Wx2j[pos + 2*tsize] = &y2j[iwest[pos + 2*tsize]];
                // Wx2j[pos + 3*tsize] = &y2j[iwest[pos + 3*tsize]];
                // Wx2j[pos + 4*tsize] = &y2j[iwest[pos + 4*tsize]];
                // Wx2j[pos + 5*tsize] = &y2j[iwest[pos + 5*tsize]];
            }

            if (j == 0)
            {
                inorth[pos] = 4*tsize + nlat*nlat - i - 1;
                inorth[pos + tsize] = nlat*nlat - i - 1;
                inorth[pos + 2*tsize] = 5*tsize + nlat*nlat - i - 1;
                inorth[pos + 3*tsize] = 2*tsize + nlat*nlat - i - 1;
                inorth[pos + 4*tsize] = 1*tsize + nlat*nlat - i - 1;
                inorth[pos + 5*tsize] = 3*tsize + nlat*nlat - i - 1;
                
                unorth[pos] = &v[inorth[inorth[pos]]];
                unorth[pos + tsize] = &v[inorth[inorth[pos + tsize]]];
                unorth[pos + 2*tsize] = &v[inorth[inorth[pos + 2*tsize]]];
                unorth[pos + 3*tsize] = &v[inorth[inorth[pos + 3*tsize]]];
                unorth[pos + 4*tsize] = &v[inorth[inorth[pos + 4*tsize]]];
                unorth[pos + 5*tsize] = &v[inorth[inorth[pos + 5*tsize]]];

                vnorth[pos] = &u[inorth[pos]];
                vnorth[pos + tsize] = &u[inorth[pos + tsize]];
                vnorth[pos + 2*tsize] = &u[inorth[pos + 2*tsize]];
                vnorth[pos + 3*tsize] = &u[inorth[pos + 3*tsize]];
                vnorth[pos + 4*tsize] = &u[inorth[pos + 4*tsize]];
                vnorth[pos + 5*tsize] = &u[inorth[pos + 5*tsize]];

                Ny1i[pos] = &x2j[inorth[pos]%tsize];
                // Ny1i[pos + tsize] = &x2j[inorth[pos + tsize]];
                // Ny1i[pos + 2*tsize] = &x2j[inorth[pos + 2*tsize]];
                // Ny1i[pos + 3*tsize] = &x2j[inorth[pos + 3*tsize]];
                // Ny1i[pos + 4*tsize] = &x2j[inorth[pos + 4*tsize]];
                // Ny1i[pos + 5*tsize] = &x2j[inorth[pos + 5*tsize]];
                Ny1j[pos] = &x2i[inorth[pos]%tsize];
                // Ny1j[pos + tsize] = &x2i[inorth[pos + tsize]];
                // Ny1j[pos + 2*tsize] = &x2i[inorth[pos + 2*tsize]];
                // Ny1j[pos + 3*tsize] = &x2i[inorth[pos + 3*tsize]];
                // Ny1j[pos + 4*tsize] = &x2i[inorth[pos + 4*tsize]];
                // Ny1j[pos + 5*tsize] = &x2i[inorth[pos + 5*tsize]];

                if (i == nlat - 1)
                {
                    unorth[pos] = &u[nlat*(nlat - 1) + tsize];
                    unorth[pos + tsize] = &u[nlat*(nlat - 1) + 4*tsize];
                    unorth[pos + 2*tsize] = &u[nlat*(nlat - 1) + 3*tsize];
                    unorth[pos + 3*tsize] = &u[nlat*(nlat - 1) + 5*tsize];
                    unorth[pos + 4*tsize] = &u[nlat*(nlat - 1)];
                    unorth[pos + 5*tsize] = &u[nlat*(nlat - 1) + 2*tsize];
                }
            }

            if (j == nlat - 1)
            {
                isouth[pos] = 5*tsize + i;
                isouth[pos + tsize] = 2*tsize + i;
                isouth[pos + 2*tsize] = 4*tsize + i;
                isouth[pos + 3*tsize] = i;
                isouth[pos + 4*tsize] = 3*tsize + i;
                isouth[pos + 5*tsize] = tsize + i;

                usouth[pos] = &v[isouth[pos]];
                usouth[pos + tsize] = &v[isouth[pos + tsize]];
                usouth[pos + 2*tsize] = &v[isouth[pos + 2*tsize]];
                usouth[pos + 3*tsize] = &v[isouth[pos + 3*tsize]];
                usouth[pos + 4*tsize] = &v[isouth[pos + 4*tsize]];
                usouth[pos + 5*tsize] = &v[isouth[pos + 5*tsize]];

                vsouth[pos] = &u[isouth[pos]];
                vsouth[pos + tsize] = &u[isouth[pos + tsize]];
                vsouth[pos + 2*tsize] = &u[isouth[pos + 2*tsize]];
                vsouth[pos + 3*tsize] = &u[isouth[pos + 3*tsize]];
                vsouth[pos + 4*tsize] = &u[isouth[pos + 4*tsize]];
                vsouth[pos + 5*tsize] = &u[isouth[pos + 5*tsize]];
            }
        }
    }

    for (int p = 0; p < 6; ++p)
    {
        int tsize = nlat*nlat;

        for (int i = 0; i < nlat; ++i)
        {
            for (int j = 0; j < nlat; ++j)
             {
                int pos = i*nlat + j;
                inortheast[pos + p*tsize] = ieast[inorth[pos + p*tsize]];
                inorthwest[pos + p*tsize] = iwest[inorth[pos + p*tsize]];
                isoutheast[pos + p*tsize] = ieast[isouth[pos + p*tsize]];
                isouthwest[pos + p*tsize] = iwest[isouth[pos + p*tsize]];
                unorthwest[pos + p*tsize] = &u[inorthwest[pos + p*tsize]];
                usouthwest[pos + p*tsize] = &u[isouthwest[pos + p*tsize]];
                vnortheast[pos + p*tsize] = &v[inortheast[pos + p*tsize]];
                vnorthwest[pos + p*tsize] = &v[inorthwest[pos + p*tsize]];

                if (j == 0)
                {
                    unorth[pos] = &v[inorth[inorth[pos]]];
                    unorth[pos + tsize] = &v[inorth[inorth[pos + tsize]]];
                    unorth[pos + 2*tsize] = &v[inorth[inorth[pos + 2*tsize]]];
                    unorth[pos + 3*tsize] = &v[inorth[inorth[pos + 3*tsize]]];
                    unorth[pos + 4*tsize] = &v[inorth[inorth[pos + 4*tsize]]];
                    unorth[pos + 5*tsize] = &v[inorth[inorth[pos + 5*tsize]]];

                    inorthwest[pos + p*tsize] = isouth[inorth[pos + p*tsize]];
                    inortheast[pos + p*tsize] = inorth[inorth[pos + p*tsize]];
                    unorthwest[pos + p*tsize] = &v[inorth[inorthwest[pos + p*tsize]]];
                    vnortheast[pos + p*tsize] = &u[inortheast[pos + p*tsize]];
                    vnorthwest[pos + p*tsize] = &u[inorthwest[pos + p*tsize]];
                }

                if (i == 0)
                {
                    inorthwest[pos + p*tsize] = iwest[iwest[pos + p*tsize]];
                    isouthwest[pos + p*tsize] = ieast[iwest[pos + p*tsize]];
                    unorthwest[pos + p*tsize] = &v[inorthwest[pos + p*tsize]];
                    usouthwest[pos + p*tsize] = &v[isouthwest[pos + p*tsize]];
                    vnorthwest[pos + p*tsize] = &u[inorthwest[pos + p*tsize]];
                }

                if (j == nlat - 1)
                {
                    isoutheast[pos + p*tsize] = isouth[isouth[pos + p*tsize]];
                    isouthwest[pos + p*tsize] = inorth[isouth[pos + p*tsize]];
                    usouthwest[pos + p*tsize] = &v[isouthwest[pos + p*tsize]];
                }

                if (i == nlat - 1)
                {
                    veast[pos] = &u[iwest[ieast[pos]]];
                    veast[pos + tsize] = &u[iwest[ieast[pos + tsize]]];
                    veast[pos + 2*tsize] = &u[iwest[ieast[pos + 2*tsize]]];
                    veast[pos + 3*tsize] = &u[iwest[ieast[pos + 3*tsize]]];
                    veast[pos + 4*tsize] = &u[iwest[ieast[pos + 4*tsize]]];
                    veast[pos + 5*tsize] = &u[iwest[ieast[pos + 5*tsize]]];

                    inortheast[pos + p*tsize] = ieast[ieast[pos + p*tsize]];
                    isoutheast[pos + p*tsize] = iwest[ieast[pos + p*tsize]];
                    vnortheast[pos + p*tsize] = &u[iwest[inortheast[pos + p*tsize]]];
                }
            }
        }

        inorthwest[p*tsize] = iwest[p*tsize];
        isouthwest[nlat - 1 + p*tsize] = isouth[nlat - 1 + p*tsize];
        inortheast[nlat*(nlat - 1) + p*tsize] = inorth[nlat*(nlat - 1) + p*tsize];
        isoutheast[nlat*nlat - 1 + p*tsize] = ieast[nlat*nlat - 1 + p*tsize];

        unorthwest[p*tsize] = &v[inorthwest[p*tsize]];
        usouthwest[nlat - 1 + p*tsize] = &v[isouthwest[nlat - 1 + p*tsize]];
        vnortheast[nlat*(nlat - 1) + p*tsize] = &u[inortheast[nlat*(nlat - 1) + p*tsize]];
        vnorthwest[p*tsize] = &u[inorthwest[p*tsize]];
    }

    // std::remove("u.csv");
    // std::remove("v.csv");
    // std::remove("s.csv");
    // std::remove("p.csv");
    // std::remove("t.csv");
    // std::remove("g.csv");

    // std::ofstream mdata("metadata.csv");
    // std::ostringstream streamObj;
    // std::string strObj;
    // streamObj << nlat;
    // strObj = streamObj.str();
    // mdata << strObj << ',';
    // streamObj.str(std::string());
    // streamObj << nlevel;
    // strObj = streamObj.str();
    // mdata << strObj;
    // mdata.close();
    // strObj.clear();

    // ps[nlat*(nlat/2) + nlat/2] += 2000;
    // ps_buffer[nlat*(nlat/2) + nlat/2] = ps[nlat*(nlat/2) + nlat/2];

    for (int i = 0; i < lsize; ++i)
    {
        double pi = ps[i] - ptop;

        for (int level = 0; level < nlevel - 1; ++level)
        {
            int l_ = 6*nlat*nlat*level;

            theta[i + l_] = Ts*pi;
            theta_buffer[i + l_] = theta[i + l_];

            double p = ((sigma[level] + ds/2)*pi + ptop);
            T[i + l_] = theta[i + l_]*pow((p_ref/p), -R/cp)/pi;
            friction[i + l_] = level < nlevel - 2 ? 0 : (1 + terrain[i]/10)*0.000015*dt;
        }
    }

    // theta[nlat*(nlat/2) + nlat/2 + lsize*(nlevel - 2)] = 320*(ps[nlat*(nlat/2) + nlat/2] - ptop);
    // theta_buffer[nlat*(nlat/2) + nlat/2 + lsize*(nlevel - 2)] = theta[nlat*(nlat/2) + nlat/2 + lsize*(nlevel - 2)];

    Hydrostatic();
}

SkyCube::~SkyCube()
{
    delete latitude;
    delete longitude;
    delete sigma;
    delete area;
    delete ps;
    delete ps_buffer;
    delete ps_dt;
    delete surface_pot;
    delete x1i;
    delete x1j;
    delete x2i;
    delete x2j;
    delete y1i;
    delete y1j;
    delete y2i;
    delete y2j;
    delete Wx2i;
    delete Wx2j;
    delete Ny1i;
    delete Ny1j;
    delete eastlength;
    delete westlength;
    delete northlength;
    delete southlength;
    delete u;
    delete u_buffer;
    delete u_accel;
    delete v;
    delete v_buffer;
    delete v_accel;
    delete theta;
    delete theta_buffer;
    delete theta_dt;
    delete T;
    delete geopot;
    delete sigma_vel;
    delete friction;
    delete inorth;
    delete isouth;
    delete ieast;
    delete iwest;
    delete inortheast;
    delete inorthwest;
    delete isoutheast;
    delete isouthwest;
    delete unorth;
    delete usouth;
    delete ueast;
    delete uwest;
    delete vnorth;
    delete vsouth;
    delete veast;
    delete vwest;
    delete unorthwest;
    delete usouthwest;
    delete vnortheast;
    delete vnorthwest;
}