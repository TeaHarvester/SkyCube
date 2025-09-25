#include "Atmosphere.h"
#include <cmath>
#include <iostream>
#include <iomanip>

void Atmosphere::Dynamics()
{
    SolvePrimitive();
    EulerForward();

    double iter = 0;

    while (iter < 1000)
    {
        // WriteState();

        for (int i = 0; i < 100; ++i)
        {
            SolvePrimitive();
            LeapFrog();
        }

        // PrintV(1);
        Print(1, 't');

        ++iter;
    }
}

void Atmosphere::SolvePrimitive()
{
    Hydrostatic();
    SurfacePressureTendency();
    SigmaVelocity();
    ThetaTendency();

    double k = R/cp;

    for (int level = 0; level < nlevel - 1; ++level)
    {
        for (int i = 0; i < xsize*ysize; ++i)
        {
            double _pi[4];
            double _phi[4];
            double _theta[4];

            _pi[0] = ps[i] - ptop;
            _pi[1] = ps[ieast[i]] - ptop;
            _pi[2] = ps[i] - ptop;
            _pi[3] = ps[isouth[i]] - ptop;
            _phi[0] = geopot[level][i];
            _phi[1] = geopot[level][ieast[i]];
            _phi[2] = geopot[level][i];
            _phi[3] = geopot[level][isouth[i]];
            _theta[0] = theta[level][i];
            _theta[1] = theta[level][ieast[i]];
            _theta[2] = theta[level][i];
            _theta[3] = theta[level][isouth[i]];

            double pix = (_pi[0] + _pi[1])/2;
            double piy = (_pi[2] + _pi[3])/2;

            std::vector<double> advect = Upstream(level, i);

            double sigma_advect[2] = {0, 0};
            double upper_u = level == 0 ? u[level][i] : u[level - 1][i];
            double upper_v = level == 0 ? v[level][i] : v[level - 1][i];
            double lower_u = level == nlevel - 2 ? u[level][i] : u[level + 1][i];
            double lower_v = level == nlevel - 2 ? v[level][i] : v[level + 1][i];

            double sigma_velx = (sigma_vel[level][i] + sigma_vel[level][ieast[i]] + 
                                sigma_vel[level + 1][i] + sigma_vel[level + 1][ieast[i]])/4;
            double sigma_vely = (sigma_vel[level][i] + sigma_vel[level][isouth[i]] + 
                                sigma_vel[level + 1][i] + sigma_vel[level + 1][isouth[i]])/4;

            double s_velx_f = (sigma_velx + fabs(sigma_velx))/2;
            double s_velx_b = (sigma_velx - fabs(sigma_velx))/2;
            double s_vely_f = (sigma_vely + fabs(sigma_vely))/2;
            double s_vely_b = (sigma_vely - fabs(sigma_vely))/2;

            sigma_advect[0] = (s_velx_f*(u[level][i] - lower_u) + s_velx_b*(upper_u - u[level][i]))/ds;
            sigma_advect[1] = (s_vely_f*(v[level][i] - lower_v) + s_vely_b*(upper_v - v[level][i]))/ds;

            double pressure_grad[2];
            double px = (sigma[level] + ds/2)*pix + ptop;
            double p1x = sigma[level]*pix + ptop;
            double p2x = sigma[level + 1]*pix + ptop;
            double py = (sigma[level] + ds/2)*piy + ptop;
            double p1y = sigma[level]*piy + ptop;
            double p2y = sigma[level + 1]*piy + ptop;
            // double Tx = (((_theta[0] + _theta[1])/2)/((k + 1)*pow(p_ref, k)*ds*pix*pix))*(pow(p2x, k + 1) - pow(p1x, k + 1));
            // double Ty = (((_theta[2] + _theta[3])/2)/((k + 1)*pow(p_ref, k)*ds*piy*piy))*(pow(p2y, k + 1) - pow(p1y, k + 1));
            double Tx = (((_theta[0] + _theta[1])/2)/((k + 1)*pow(p_ref, k)*ds*pix))*(pow(p2x, k + 1) - pow(p1x, k + 1));
            double Ty = (((_theta[2] + _theta[3])/2)/((k + 1)*pow(p_ref, k)*ds*piy))*(pow(p2y, k + 1) - pow(p1y, k + 1));
            
            double rhox = px/(R*Tx);
            double rhoy = py/(R*Ty);
            pressure_grad[0] = (pix/rhox)*(sigma[level] + ds/2)*(_pi[1] - _pi[0])/dx[i/ysize];
            pressure_grad[1] = (piy/rhoy)*(sigma[level] + ds/2)*(_pi[2] - _pi[3])/dy[i%ysize];;

            double geograd[2];
            geograd[0] = pix*(_phi[1] - _phi[0])/dx[i/ysize];
            geograd[1] = piy*(_phi[2] - _phi[3])/dy[i%ysize];

            u_accel[level][i] = -advect[0] - sigma_advect[0] - geograd[0] - pressure_grad[0];
            v_accel[level][i] = -advect[1] - sigma_advect[1] - geograd[1] - pressure_grad[1];
        }
    }
}

std::vector<double> Atmosphere::Upstream(int& level, int& pos)
{
    int& i = pos;
    double _pi[4];
    double _u[4];
    double _v[4];

    _pi[0] = ps[i] - ptop;
    _pi[1] = ps[ieast[i]] - ptop;
    _pi[2] = ps[i] - ptop;
    _pi[3] = ps[isouth[i]] - ptop;
    _u[0] = u[level][iwest[i]];
    _u[1] = u[level][ieast[i]];
    _u[2] = u[level][inorth[i]];
    _u[3] = u[level][isouth[i]];
    _v[0] = v[level][iwest[i]];
    _v[1] = v[level][ieast[i]];
    _v[2] = v[level][inorth[i]]*ipolar[inorth[i]];
    _v[3] = v[level][isouth[i]]*ipolar[isouth[i]];

    double upper_v = (_v[2] + v[level][inortheast[i]]*ipolar[inortheast[i]])/2;
    double lower_v = (v[level][i] + _v[1])/2;
    double west_u = (_u[0] + u[level][isouthwest[i]])/2;
    double east_u = (u[level][i] + _u[3])/2;

    double upper_pi = (ps[ieast[i]] + ps[i] + ps[inortheast[i]] + ps[inorth[i]] - 4*ptop)/4;
    double lower_pi = (ps[ieast[i]] + ps[i] + ps[isoutheast[i]] + ps[isouth[i]] - 4*ptop)/4;
    double west_pi = (ps[iwest[i]] + ps[i] + ps[isouthwest[i]] + ps[isouth[i]] - 4*ptop)/4;
    double east_pi = (ps[ieast[i]] + ps[i] + ps[isoutheast[i]] + ps[isouth[i]] - 4*ptop)/4;

    double ux = (u[level][i] + (_u[0]/2 + _u[1]/2))/2;
    double vy = (v[level][i] + (_v[2]/2 + _v[3]/2))/2;

    double ux_f = (ux + fabs(ux))/(_pi[0] + _pi[1]);
    double ux_b = (ux - fabs(ux))/(_pi[0] + _pi[1]);
    double vx_f = (upper_v + lower_v + fabs(upper_v + lower_v))/(2*(upper_pi + lower_pi));
    double vx_b = (upper_v + lower_v - fabs(upper_v + lower_v))/(2*(upper_pi + lower_pi));

    double vy_f = (vy + fabs(vy))/(_pi[2] + _pi[3]);
    double vy_b = (vy - fabs(vy))/(_pi[2] + _pi[3]);
    double uy_f = (west_u + east_u + fabs(west_u + east_u))/(2*(west_pi + east_pi));
    double uy_b = (west_u + east_u - fabs(west_u + east_u))/(2*(west_pi + east_pi));

    double u_advect = (ux_f*(u[level][i] - _u[0]) + ux_b*(_u[1] - u[level][i]))/dx[i/ysize] + 
                      (vx_f*(u[level][i] - _u[3]) + vx_b*(_u[2] - u[level][i]))/dy[i%ysize];

    double v_advect = (vy_f*(v[level][i] - _v[3]) + vy_b*(_v[2] - v[level][i]))/dy[i%ysize] + 
                      (uy_f*(v[level][i] - _v[0]) + uy_b*(_v[1] - v[level][i]))/dx[i/ysize];

    std::vector<double> advect = {u_advect, v_advect};
    return advect;
}

void Atmosphere::Hydrostatic()
{
    double k = R/cp;

    for (int i = 0; i < xsize*ysize; ++i)
    { 
        double pi = ps[i] - ptop;
        double potential = 0;

        for (int level = 0; level < nlevel - 2; ++level)
        {
            double p0 = pi*sigma[level] + ptop;
            double p = pi*sigma[level + 1] + ptop;
            double p1 = pi*sigma[level + 2] + ptop;
            double PI1 = (cp/((k + 1)*pow(p_ref, k)*pi*ds))*(pow(p, k + 1) - pow(p0, k + 1));
            double PI2 = (cp/((k + 1)*pow(p_ref, k)*pi*ds))*(pow(p1, k + 1) - pow(p, k + 1));
            double PI_dpi = ((cp/(pow(p_ref, k)*ds))*(sigma[level + 1]*pow(p, k) - sigma[level]*pow(p0, k)) - PI1)/pi;
            double theta_inter = (theta[level][i] + theta[level + 1][i])/2;

            potential += pi*theta[level][i]*PI_dpi*ds;
            potential -= sigma[level + 1]*(PI2 - PI1)*theta_inter;
        }

        int level = nlevel - 2;
        double p0 = pi*sigma[level] + ptop;
        double p = pi*sigma[level + 1] + ptop;
        double PI1 = (cp/((k + 1)*pow(p_ref, k)*pi*ds))*(pow(p, k + 1) - pow(p0, k + 1));
        double PI_dpi = ((cp/(pow(p_ref, k)*ds))*(sigma[level + 1]*pow(p, k) - sigma[level]*pow(p0, k)) - PI1)/pi;

        potential += pi*theta[level][i]*PI_dpi*ds;
        geopot[level][i] = surface_pot[i] + potential;

        for (int level = nlevel - 3; level >= 0; --level)
        {
            double p0 = pi*sigma[level] + ptop;
            double p = pi*sigma[level + 1] + ptop;
            double p1 = pi*sigma[level + 2] + ptop;
            double PI1 = (cp/((k + 1)*pow(p_ref, k)*pi*ds))*(pow(p, k + 1) - pow(p0, k + 1));
            double PI2 = (cp/((k + 1)*pow(p_ref, k)*pi*ds))*(pow(p1, k + 1) - pow(p, k + 1));
            double theta_inter = (theta[level][i] + theta[level + 1][i])/2;
            geopot[level][i] = geopot[level + 1][i] + theta_inter*(PI2 - PI1);
        }
    }  
}

void Atmosphere::SurfacePressureTendency()
{
    double _u[4];
    double _v[4];

    for (int i = 0; i < xsize*ysize; ++i)
    {
        double ps_vel = 0;

        for(int level = 0; level < nlevel - 1; ++level)
        {
            _u[0] = u[level][iwest[i]];
            _u[1] = u[level][i];
            _v[0] = v[level][inorth[i]]*ipolar[inorth[i]];
            _v[1] = v[level][i];

            ps_vel -= ds*(_u[1] - _u[0])/dx[i/ysize];
            ps_vel -= ds*(_v[0] - _v[1])/dy[i%ysize];
        }

        ps_dt[i] = ps_vel;
    }
}

void Atmosphere::SigmaVelocity()
{
    double _u[2];
    double _v[2];

    for (int i = 0; i < xsize*ysize; ++i)
    {
        double ps_vel = 0;

        for (int level = 0; level < nlevel - 1; ++level)
        {
            double pi = ps[i] - ptop;
            _u[0] = u[level][iwest[i]];
            _u[1] = u[level][i];
            _v[0] = v[level][inorth[i]]*ipolar[inorth[i]];
            _v[1] = v[level][i];

            ps_vel -= ds*(_u[1] - _u[0])/dx[i/ysize];
            ps_vel -= ds*(_v[0] - _v[1])/dy[i%ysize];
            sigma_vel[level + 1][i] = (sigma[level + 1]*ps_dt[i] - ps_vel)/pi;
        }
    }
} 

void Atmosphere::ThetaTendency()
{
    // Radiate();
    double _theta[4];
    double _u[2];
    double _v[2];

    for (int i = 0; i < xsize*ysize; ++i)
    {
        for (int level = 0; level < nlevel - 1; ++level)
        {
            double pi = ps[i] - ptop;
            _theta[0] = theta[level][iwest[i]];
            _theta[1] = theta[level][ieast[i]];
            _theta[2] = theta[level][inorth[i]];
            _theta[3] = theta[level][isouth[i]];
            _u[0] = u[level][iwest[i]];
            _u[1] = u[level][i];
            _v[0] = v[level][inorth[i]]*ipolar[inorth[i]];
            _v[1] = v[level][i];
            // _pi[0] = ps[dir[0]] - ptop;
            // _pi[1] = ps[dir[1]] - ptop;
            // _pi[2] = ps[dir[2]] - ptop;
            // _pi[3] = ps[dir[4]] - ptop;


            double ux_f = (_u[0] + _u[1] + fabs(_u[0] + _u[1]))/4;
            double ux_b = (_u[0] + _u[1] - fabs(_u[0] + _u[1]))/4;
            double vy_f = (_v[0] + _v[1] + fabs(_v[0] + _v[1]))/4;
            double vy_b = (_v[0] + _v[1] - fabs(_v[0] + _v[1]))/4;


            double advect = (ux_f*(theta[level][i] - _theta[0]) + 
                             ux_b*(_theta[1] - theta[level][i]))/dx[i/ysize] + 
                            (vy_f*(theta[level][i] - _theta[3]) +
                             vy_b*(_theta[2] - theta[level][i]))/dy[i%ysize];

            double sigma_vel_avg = (sigma_vel[level][i] + sigma_vel[level + 1][i])/2;
            double s_vel_f = (sigma_vel_avg + fabs(sigma_vel_avg))/2;
            double s_vel_b = (sigma_vel_avg - fabs(sigma_vel_avg))/2;

            double upper_theta = level == 0 ? theta[level][i] : theta[level - 1][i];
            double lower_theta = level == nlevel - 2 ? theta[level][i] : theta[level + 1][i];

            double sigma_advect = pi*(s_vel_f*(theta[level][i] - lower_theta) + s_vel_b*(upper_theta - theta[level][i]))/ds;
            double heatflux = pi*(theta[level][i]/T[level][i])*(radiant_heatflux[level][i])/cp;
            theta_dt[level][i] = (heatflux - advect - sigma_advect)/pi;
        }
    }
}

void Atmosphere::Radiate()
{
    double irradiance = 340;
    double em = 0.78;

    // std::random_device rd; 
    // std::mt19937 gen(rd());
    // std::uniform_int_distribution<> distrib(-10, 10);

    for (int i = 0; i < xsize*ysize; ++i)
    {
        // double incoming = irradiance + (double)distrib(gen);
        double incoming = irradiance*(1 - albedo[i])*cos(latitude[i%ysize]);
        double outgoing = (1 - em)*boltz*pow(T[nlevel - 2][i], 4);
        double rho = ((1 - ds/2)*(ps[i] - ptop) + ptop)/(R*T[nlevel - 2][i]);
        double dphi = (ps[i] - ptop)*ds/rho;
        radiant_heatflux[nlevel - 2][i] = (incoming - outgoing)*g/(rho*dphi);

        for (int leve = 0; leve < nlevel - 2; ++leve)
        {
            outgoing = em*boltz*pow(T[leve][i], 4);
            rho = ((1 - ds/2)*(ps[i] - ptop) + ptop)/(R*T[leve][i]);
            dphi = (ps[i] - ptop)*ds/rho;
            radiant_heatflux[leve][i] = -outgoing*g/(rho*dphi);
        }
    }
}

void Atmosphere::LeapFrog()
{
    double p_swap;
    double u_swap;
    double v_swap;
    double theta_swap;
    double timefilter = 0.01;

    for (int i = 0; i < xsize*ysize; ++i)
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
        for (int i = 0; i < xsize*ysize; ++i)
        {
            u_swap = u[level][i];
            v_swap = v[level][i];
            theta_swap = theta[level][i];
            u[level][i] = u_buffer[level][i] + u_accel[level][i]*dt;
            v[level][i] = v_buffer[level][i] + v_accel[level][i]*dt;
            theta[level][i] = theta_buffer[level][i] + theta_dt[level][i]*dt;
            T[level][i] = theta[level][i]*pow((p_ref/((sigma[level] + ds/2)*(ps[i] - ptop) + ptop)), -R/cp);

            // implicit schemes
            double f = 2*OM*std::sin(latitude[i%ysize]);

            double u_interp = (u[level][i] + u[level][iwest[i]]
             + u[level][isouth[i]] + u[level][isouthwest[i]])/4;

            double v_interp = (v[level][i] + v[level][inorth[i]]
             + v[level][ieast[i]] + v[level][inortheast[i]])/4;

            u[level][i] = (u[level][i] - f*dt*v_interp)/(pow(f*dt, 2) + friction[level][i] + 1);
            v[level][i] = (v[level][i] + f*dt*u_interp)/(pow(f*dt, 2) + friction[level][i] + 1);

            u_buffer[level][i] = u_swap + timefilter*(u[level][i] -2*u_swap + u_buffer[level][i]);
            v_buffer[level][i] = v_swap + timefilter*(v[level][i] -2*v_swap + v_buffer[level][i]);
            theta_buffer[level][i] = theta_swap + timefilter*(theta[level][i] -2*theta_swap + theta_buffer[level][i]);
        }
    }
}

void Atmosphere::EulerForward()
{
    for (int i = 0; i < xsize*ysize; ++i)
    {
        ps[i] = ps[i] + ps_dt[i]*(dt/2); 
    }

    for (int level = 0; level < nlevel - 1; ++level)
    {
        for (int i = 0; i < xsize*ysize; ++i)
        {
            u[level][i] = u[level][i] + u_accel[level][i]*(dt/2);
            v[level][i] = v[level][i] + v_accel[level][i]*(dt/2);
        }   
    }

    ThetaTendency();

    for (int level = 0; level < nlevel - 1; ++level)
    {
        for (int i = 0; i < xsize*ysize; ++i)
        {
            theta[level][i] = theta[level][i] + theta_dt[level][i]*(dt/2);
            T[level][i] = theta[level][i]*pow((p_ref/((sigma[level] + ds/2)*(ps[i] - ptop) + ptop)), -R/cp);
        }
    }
}

std::vector<int> Atmosphere::GetPeriodicDir(const int centre)
{
    const int xcentre = centre/ysize;
    const int ycentre = centre%ysize;

    std::vector<int> dir(6);
    dir[0] = xcentre > 0 ? ycentre + ysize*((xcentre - 1)) : ycentre + ysize*(xsize - 1);
    dir[1] = xcentre < xsize - 1 ? ycentre + ysize*((xcentre + 1)) : ycentre;
    dir[2] = ycentre > 0 ? centre - 1 : xcentre*ysize + (ysize - 1);
    // dir[3] = dir[2] + ysize*st;
    dir[4] = ycentre < ysize - 1 ? centre + 1 : xcentre*ysize;
    // dir[5] = dir[4] + ysize*st;

    return dir;
}

std::vector<int> Atmosphere::GetMercatorDir(const int centre)
{
    const int xcentre = centre/ysize;
    const int ycentre = centre%ysize;

    std::vector<int> dir(6);
    dir[0] = xcentre > 0 ? ycentre + ysize*((xcentre - 1)) : ycentre + ysize*(xsize - 1);
    dir[1] = xcentre < xsize - 1 ? ycentre + ysize*((xcentre + 1)) : ycentre;
    dir[2] = ycentre > 0 ? centre - 1 : ysize*((xcentre + xsize/2)%xsize);
    // dir[3] = dir[2] + ysize*st;
    dir[4] = ycentre < ysize - 1 ? centre + 1 : ysize*((xcentre + xsize/2)%xsize) + ysize - 1;
    // dir[5] = dir[4] + ysize*st;

    return dir;
}

void Atmosphere::PrintV(const int level)
{
    std::cout << std::fixed;
    std::cout << std::endl;

    for (int i = 0; i < 2*ysize; ++i)
    {
        std::cout << std::endl;

        for (int j = 0; j < xsize; ++j)
        {
            if (i%2 == 0)
            {
                std::cout << std::setprecision(0) << ps[i/2 + j*ysize]/100 << " " << 
                // std::cout << std::setprecision(0) << geopot[level][i/2 + j*size] << " " << 
                // std::cout << std::setprecision(2) << sigma_vel[level][i/2 + j*size] << ' ' << 
                std::setprecision(2) << u[level][i/2 + j*ysize] << " ";
            }

            else
            {
                std::cout << v[level][i/2 + j*ysize] << "      ";
            }
        }
    }
}

void Atmosphere::Print(const int level, char var)
{
    std::vector<std::vector<double>>* output;
    std::vector<double> p(xsize*ysize, 2);
    double q = 1;

    switch (var)
    {
    case 'u':
        output = &u;
        p = ps;
        q = ptop;
        break;
    case 'v':
        output = &v;
        p = ps;
        q = ptop;
        break;
    case 'g':
        output = &geopot;
        break;
    case 'T':
        output = &T;
        break;
    case 't':
        output = &theta;
        // p = ps;
        // q = ptop;
        break;
    case 's':
        output = &sigma_vel;
        break;
    default:
        return;
    }

    std::cout << std::fixed << std::setprecision(2);
    std::cout << std::endl;

    for (int i = 0; i < ysize; ++i)
    {
        std::cout << std::endl;

        for (int j = 0; j < xsize; ++j)
        {
            std::cout << (*output)[level][i + j*ysize]/(p[i + j*ysize] - q) << ' ';
        }
    }
}

void Atmosphere::PrintP(char var)
{
    std::vector<double>* output;

    switch (var)
    {
    case 's':
        output = &ps;
        break;
    case 'v':
        output = &ps_dt;
        break;
    default:
        return;
    }


    std::cout << std::fixed << std::setprecision(2);
    std::cout << std::endl;

    for (int i = 0; i < ysize; ++i)
    {
        std::cout << std::endl;

        for (int j = 0; j < xsize; ++j)
        {
            std::cout << (*output)[i + j*ysize]/100 << ' ';
        }
    }
}

void Atmosphere::WriteState()
{
    std::ofstream ufile("u.csv", std::ios_base::app);
    std::ofstream vfile("v.csv", std::ios_base::app);
    std::ofstream sfile("s.csv", std::ios_base::app);
    std::ofstream pfile("p.csv", std::ios_base::app);
    std::ofstream gfile("g.csv", std::ios_base::app);
    std::ofstream tfile("t.csv", std::ios_base::app);
    std::ostringstream streamObj;
    std::string strObj;
    

    for (int level = 0; level < nlevel - 1; ++level)
    {
        for (int i = 0; i < xsize*ysize; ++i)
        {
            streamObj << std::scientific << std::setprecision(3);
            streamObj << (u[level][i] + u[level][iwest[i]])/(2*(ps[i] - ptop));
            strObj = streamObj.str();
            ufile << strObj << ',';
            streamObj.str(std::string());

            streamObj << (v[level][i] + v[level][inorth[i]])/(2*(ps[i] - ptop));
            strObj = streamObj.str();
            vfile << strObj << ',';
            streamObj.str(std::string());

            streamObj << sigma_vel[level][i];
            strObj = streamObj.str();
            sfile << strObj << ',';
            streamObj.str(std::string());

            streamObj << geopot[level][i];
            strObj = streamObj.str();
            gfile << strObj << ',';
            streamObj.str(std::string());

            // streamObj << theta[level][i]/(ps[i] - ptop);
            streamObj << std::scientific << std::setprecision(7);
            streamObj << theta[level][i];
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

    streamObj << std::scientific << std::setprecision(5);

    for (int i = 0; i < xsize*ysize; ++i)
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

Atmosphere::Atmosphere(int nlat, int n_lev, double R, double p0, double pt, std::vector<int>& terre)
:
xsize(6*nlat),
ysize(nlat),
nlevel(n_lev),
dt(0.001),
ptop(pt),
p_ref(p0),
OM(0.000072921),
radius(R),
ds(1/((double)n_lev - 1)),
terrain(terre)
{
}

Atmosphere::Atmosphere(int sx, int sy, int n_lev, double dy_, double p0, double pt, double lt, std::vector<int>& terre)
:
xsize(sx),
ysize(sy),
nlevel(n_lev),
dt(0.01),
ptop(pt),
p_ref(p0),
OM(0.000072921),
radius(6378000),
ds(1/((double)n_lev - 1)),
dx(xsize, dy_),
dy(ysize, dy_),
ps(sx*sy),
ps_buffer(sx*sy),
ps_dt(sx*sy, 0),
surface_pot(sx*sy, 0),
albedo(sx*sy, 0),
latitude(sy, lt),
T(n_lev - 1, std::vector<double>(sx*sy)),
theta(n_lev - 1, std::vector<double>(sx*sy)),
theta_buffer(n_lev - 1, std::vector<double>(sx*sy)),
theta_dt(n_lev - 1, std::vector<double>(sx*sy, 0)),
geopot(n_lev - 1, std::vector<double>(sx*sy, 0)),
u(n_lev - 1, std::vector<double>(sx*sy, 0)),
v(n_lev - 1, std::vector<double>(sx*sy, 0)),
u_buffer(n_lev - 1, std::vector<double>(sx*sy, 0)),
v_buffer(n_lev - 1, std::vector<double>(sx*sy, 0)),
sigma_vel(n_lev, std::vector<double>(sx*sy, 0)),
u_accel(n_lev - 1, std::vector<double>(sx*sy, 0)),
v_accel(n_lev - 1, std::vector<double>(sx*sy, 0)),
friction(n_lev - 1, std::vector<double>(sx*sy, 0)),
radiant_heatflux(n_lev - 1, std::vector<double>(sx*sy, 0)),
sigma(n_lev),
terrain(terre),
inorth(sx*sy),
isouth(sx*sy),
ieast(sx*sy),
iwest(sx*sy),
inortheast(sx*sy),
inorthwest(sx*sy),
isoutheast(sx*sy),
isouthwest(sx*sy),
ipolar(sx*sy, 1)
{
    double Ts = 298;

    for (int i = 0; i < xsize*ysize; ++i)
    {   
        surface_pot[i] = 10*terrain[i]*g;
        ps[i] = pow(pow(p_ref, R/cp)*(1 - surface_pot[i]/(cp*Ts)), cp/R);
        ps_buffer[i] = ps[i];
        albedo[i] = terrain[i] > 24 ? 0.7 : 0;
        albedo[i] = terrain[i] == 0 ? 0.05 : 0;

        std::vector<int> dir;
        std::vector<int> dirN;
        std::vector<int> dirS;
        
        dir = GetPeriodicDir(i);
        dirN = GetPeriodicDir(dir[2]);
        dirS = GetPeriodicDir(dir[4]);
        
        inorth[i] = dir[2];
        isouth[i] = dir[4];
        ieast[i] = dir[1];
        iwest[i] = dir[0];
        inortheast[i] = dirN[1];
        inorthwest[i] = dirN[0];
        isoutheast[i] = dirS[1];
        isouthwest[i] = dirS[0];
    }

    for (int level = 0; level < nlevel; ++level)
    {
        sigma[level] = (double)level/((double)nlevel - 1);
    }

    for (int i = 0; i < xsize*ysize; ++i)
    {
        double pi = ps[i] - ptop;

        for (int level = 0; level < nlevel - 1; ++level)
        {
            theta[level][i] = Ts;
            theta_buffer[level][i] = theta[level][i];

            double p = ((sigma[level] + ds/2)*pi + ptop);
            T[level][i] = theta[level][i]*pow((p_ref/p), -R/cp);
            friction[level][i] = level > nlevel - 1 ? 0.0015*dt : (1 + terrain[i]/10)*0.0015*dt;
        }
    }

    theta[1][ysize/2 + ysize*(xsize/2)] += 0.1;
    theta_buffer[1][ysize/2 + ysize*(xsize/2)] += 0.1;
    // ps[ysize/2 + ysize*(xsize/2)] -= 1000;
    // ps[ysize/2 + ysize*(xsize/2)] -= 1000;

    Hydrostatic();
}

Atmosphere::~Atmosphere()
{}