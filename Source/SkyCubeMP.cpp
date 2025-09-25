#include "SkyCubeMP.h"
#include <iostream>
#include <iomanip>
#include <cmath>

void SkyCubeMP::Dynamics(double dt, int ntime, int write_interval)
{
    #pragma omp parallel
    {
    for (int t = 0; t < ntime; ++t)
    {
        int tsize = nlat*nlat;
        int lsize = 6*tsize;

        #pragma omp for
        for (int i = 0; i < lsize; ++i)
        { 
            // Surface pressure tendency
            int pos = i%tsize;
            int north_halo = 1 - 2*(pos%nlat == 0);
            int west_halo = 1 - 2*(pos/nlat == 0);
            double pi = ps[i] - ptop;
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

                ps_vel[level] = -ds*(ue*eastlength[pos] - uw*westlength[pos])/area[pos];
                ps_vel[level] -= ds*(vn*northlength[pos] - vs*southlength[pos])/area[pos];
                ps_dt[i] += ps_vel[level];
                cumul_ps_vel[level] = ps_dt[i]; 
                mass_change += ps_vel[level]*area[pos];

                // Horizontal heat flux
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

                // theta_dt[i + l_] = 0;

                // double _q[4];
                // _q[0] = q[iwest[i] + l_]/(ps[iwest[i]] - ptop);
                // _q[1] = q[ieast[i] + l_]/(ps[ieast[i]] - ptop);
                // _q[2] = q[inorth[i] + l_]/(ps[inorth[i]] - ptop);
                // _q[3] = q[isouth[i] + l_]/(ps[isouth[i]] - ptop);

                // // Horizontal moisture flux
                // q_dt[i + l_] = (westlength[pos]*(uw_f*_q[0] + uw_b*q[i + l_]/pi) +
                //                 southlength[pos]*(vs_f*_q[3] + vs_b*q[i + l_]/pi) -
                //                 northlength[pos]*(vn_b*_q[2] + vn_f*q[i + l_]/pi) -  
                //                 eastlength[pos]*(ue_b*_q[1] + ue_f*q[i + l_]/pi))/area[pos];
            }

            for(int level = 0; level < nlevel - 1; ++level)
            {
                int l_ = level*lsize;
                sigma_vel[i + l_ + lsize] = (sigma[level + 1]*ps_dt[i] - cumul_ps_vel[level])/pi;
            }

            double k = R/cp;
            double potential = 0;

            for (int level = 0; level < nlevel - 2; ++level)
            {
                // Hydrostatic balance
                int l_ = level*lsize;
                double p0 = pi*sigma[level] + ptop;
                double p = pi*sigma[level + 1] + ptop;
                double p1 = pi*sigma[level + 2] + ptop;
                double PI1 = (cp/((k + 1)*pow(p_ref, k)*pi*ds))*(pow(p, k + 1) - pow(p0, k + 1));
                double PI2 = (cp/((k + 1)*pow(p_ref, k)*pi*ds))*(pow(p1, k + 1) - pow(p, k + 1));
                // double PI_dpi = ((cp/(pow(p_ref, k)*ds))*(sigma[level + 1]*pow(p, k) - sigma[level]*pow(p0, k)) - PI1)/pi;
                double PI_dpi = (sigma[level] + ds/2)*cp*(pow(p, k) - pow(p0, k))/(pow(p_ref, k)*pi*ds);
                double theta_inter = (theta[i + l_] + theta[i + l_ + lsize])/2;

                potential += theta[i + l_]*pi*ds*PI_dpi;
                potential -= (sigma[level + 1])*(PI2 - PI1)*theta_inter;

                // vertical heat flux
                double& s_vel = sigma_vel[i + l_ + lsize];
                double& upper_theta = theta[i + l_];
                double& lower_theta = theta[i + l_ + lsize];
                double s_vel_f = (s_vel + fabs(s_vel))/2;
                double s_vel_b = (s_vel - fabs(s_vel))/2;
                double sigma_advect_upper = (s_vel_f*lower_theta + s_vel_b*upper_theta)/ds;
                theta_dt[i + l_] += sigma_advect_upper;
                theta_dt[i + l_ + lsize] -= sigma_advect_upper;

                // double& upper_q = q[i + l_];
                // double& lower_q = q[i + l_ + lsize];
                // sigma_advect_upper = (s_vel_f*lower_q + s_vel_b*upper_q)/ds;
                // q_dt[i + l_] += sigma_advect_upper;
                // q_dt[i + l_ + lsize] -= sigma_advect_upper;
            }

            int level = nlevel - 2;
            int l_ = level*lsize;

            double p0 = pi*sigma[level] + ptop;
            double p = pi*sigma[level + 1] + ptop;
            // double PI1 = (cp/((k + 1)*pow(p_ref, k)*pi*ds))*(pow(p, k + 1) - pow(p0, k + 1));
            double PI_dpi = (sigma[level] + ds/2)*cp*(pow(p, k) - pow(p0, k))/(pow(p_ref, k)*pi*ds);

            potential += theta[i + l_]*PI_dpi*pi*ds;
            geopot[i + l_] = surface_pot[i] + potential;

            for (int level = nlevel - 3; level >= 0; --level)
            {
                l_ = level*lsize;
                double p0 = pi*sigma[level] + ptop;
                double p = pi*sigma[level + 1] + ptop;
                double p1 = pi*sigma[level + 2] + ptop;
                double PI1 = (cp/((k + 1)*pow(p_ref, k)*pi*ds))*(pow(p, k + 1) - pow(p0, k + 1));
                double PI2 = (cp/((k + 1)*pow(p_ref, k)*pi*ds))*(pow(p1, k + 1) - pow(p, k + 1));
                double theta_inter = (theta[i + l_] + theta[i + l_ + lsize])/2;
                geopot[i + l_] = geopot[i + l_ + lsize] + theta_inter*(PI2 - PI1);
            }
        }  

        // Solve primitive equations
        #pragma omp for
        for (int i = 0; i < lsize; ++i)
        {    
            for (int level = 0; level < nlevel - 1; ++level)
            {
                int l_ = level*lsize;
                double k = R/cp;
                int pos = i%tsize;
                int epos = ieast[i]%tsize;
                int spos = isouth[i]%tsize;
                double _pi[4];
                double _PI[4];
                double _phi[4];
                double _theta[4];

                _pi[0] = ps[i] - ptop;
                _pi[1] = ps[ieast[i]] - ptop;
                _pi[2] = ps[i] - ptop;
                _pi[3] = ps[isouth[i]] - ptop;
                _PI[0] = (cp/((k + 1)*pow(p_ref, k)*_pi[0]*ds))*(pow(sigma[level + 1]*_pi[0] + ptop, k + 1) - pow(sigma[level]*_pi[0] + ptop, k + 1));
                _PI[1] = (cp/((k + 1)*pow(p_ref, k)*_pi[1]*ds))*(pow(sigma[level + 1]*_pi[1] + ptop, k + 1) - pow(sigma[level]*_pi[1] + ptop, k + 1));
                _PI[2] = (cp/((k + 1)*pow(p_ref, k)*_pi[2]*ds))*(pow(sigma[level + 1]*_pi[2] + ptop, k + 1) - pow(sigma[level]*_pi[2] + ptop, k + 1));
                _PI[3] = (cp/((k + 1)*pow(p_ref, k)*_pi[3]*ds))*(pow(sigma[level + 1]*_pi[3] + ptop, k + 1) - pow(sigma[level]*_pi[3] + ptop, k + 1));
                _phi[0] = geopot[i + l_];
                _phi[1] = geopot[ieast[i] + l_];
                _phi[2] = geopot[i + l_];
                _phi[3] = geopot[isouth[i] + l_];
                _theta[0] = theta[i + l_];
                _theta[1] = theta[ieast[i] + l_];
                _theta[2] = theta[i + l_];
                _theta[3] = theta[isouth[i] + l_];

                double pix = (_pi[0] + _pi[1])/2;
                double piy = (_pi[2] + _pi[3])/2;

                // Upstream advection
                int north_halo = 1 - 2*(pos%nlat == 0);
                int south_halo = 1 - 2*(pos%nlat == nlat - 1);
                int east_halo = 1 - 2*(pos/nlat == nlat - 1);
                int west_halo = 1 - 2*(pos/nlat == 0);
                int northeast_halo = (pos != nlat*(nlat - 1))*north_halo;
                int southwest_halo = (pos != nlat - 1)*south_halo;

                double uw = west_halo*(*(uwest[i] + l_));
                double ue = east_halo*(*(ueast[i] + l_));
                double un = *(unorth[i] + l_);
                double us = south_halo*(*(usouth[i] + l_));
                double vw = west_halo*(*(vwest[i] + l_));
                double ve = *(veast[i] + l_);
                double vn = north_halo*(*(vnorth[i] + l_));
                double vs = south_halo*(*(vsouth[i] + l_));

                double vn_interp = (vn + northeast_halo*(*(vnortheast[i] + l_)))/2;
                double vs_interp = (v[i + l_] + ve)/2;
                double west_u = (uw + southwest_halo*(*(usouthwest[i] + l_)))/2;
                double east_u = (u[i + l_] + us)/2;

                double upper_pi = (ps[ieast[i]] + ps[i] + ps[inortheast[i]] + ps[inorth[i]] - 4*ptop)/4;
                double lower_pi = (ps[ieast[i]] + ps[i] + ps[isoutheast[i]] + ps[isouth[i]] - 4*ptop)/4;
                double west_pi = (ps[iwest[i]] + ps[i] + ps[isouthwest[i]] + ps[isouth[i]] - 4*ptop)/4;
                double east_pi = (ps[ieast[i]] + ps[i] + ps[isoutheast[i]] + ps[isouth[i]] - 4*ptop)/4;

                double ux = (2*u[i + l_] + ue + uw)/2;
                double vy = (2*v[i + l_] + vn + vs)/2;

                double ux_f = (ux + fabs(ux))/(_pi[0] + _pi[1]);
                double ux_b = (ux - fabs(ux))/(_pi[0] + _pi[1]);
                double vx_f = (vn_interp + vs_interp + fabs(vn_interp + vs_interp))/(2*(upper_pi + lower_pi));
                double vx_b = (vn_interp + vs_interp - fabs(vn_interp + vs_interp))/(2*(upper_pi + lower_pi));

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

                double advect[2] = {u_advect, v_advect};
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
                // double Tx = (((_[0] + _theta[1])/2)/((k + 1)*pow(p_ref, k)*ds*pix))*(pow(p2x, k + 1) - pow(p1x, k + 1));
                // double Ty = (((_theta[2] + _theta[3])/2)/((k + 1)*pow(p_ref, k)*ds*piy))*(pow(p2y, k + 1) - pow(p1y, k + 1));
                // double rhox = px/(R*Tx);
                // double rhoy = py/(R*Ty);

                // double& u1i = x1i[pos];
                // double& u1j = x1j[pos];
                // double& u2i = x2i[pos];
                // double& u2j = x2j[pos];
                // double& v1i = y1i[pos];
                // double& v1j = y1j[pos];
                // double& v2i = y2i[pos];
                // double& v2j = y2j[pos];

                // double usin = u1i*u2j - u1j*u2i;
                // double vsin = v1i*v2j - v1j*v2i;
                // double dpi = u1i*(_pi[1] - _pi[0]);
                // double dpj = v2j*(_pi[2] - _pi[3]);
                // pressure_grad[0] = (pix/rhox)*(sigma[level] + ds/2)*dpi*usin*(2*eastlength[pos])/(area[pos] + area[epos]);
                // pressure_grad[1] = (piy/rhoy)*(sigma[level] + ds/2)*dpj*vsin*(2*southlength[pos])/(area[pos] + area[spos]);
                pressure_grad[0] = (pix*(_theta[0] + _theta[1])/2)*(_PI[1] - _PI[0])*(2*eastlength[pos])/(area[pos] + area[epos]);
                pressure_grad[1] = (piy*(_theta[2] + _theta[3])/2)*(_PI[2] - _PI[3])*(2*southlength[pos])/(area[pos] + area[spos]);

                double geograd[2];
                // double dphii = u1i*(_phi[1] - _phi[0]);
                // double dphij = v2j*(_phi[2] - _phi[3]);            
                // geograd[0] = pix*dphii*usin*(2*eastlength[pos])/(area[pos] + area[epos]);
                // geograd[1] = piy*dphij*vsin*(2*southlength[pos])/(area[pos] + area[spos]);
                geograd[0] = pix*(_phi[1] - _phi[0])*(2*eastlength[pos])/(area[pos] + area[epos]);
                geograd[1] = piy*(_phi[2] - _phi[3])*(2*southlength[pos])/(area[pos] + area[spos]);

                u_accel[i + l_] = -advect[0] - sigma_advect[0] - geograd[0] - pressure_grad[0];
                v_accel[i + l_] = -advect[1] - sigma_advect[1] - geograd[1] - pressure_grad[1];
            }
        }

        if (write_interval == -1)
        {
            EulerForward(dt);
        }

        else
        {
            LeapFrog(dt);
        }

        #pragma omp master
        {
            time += dt;
            hour_angle -= OM*dt;
            declination = tilt*sin(time/(60*60*24*365));

            if (write_interval > 0)
            {
                if (!(t%write_interval))
                {
                    WriteState();   
                }
            }
        }
    }
    }
}

void SkyCubeMP::EulerForward(double dt)
{
    int tsize = nlat*nlat;
    int lsize = 6*tsize;

    #pragma omp for
    for (int i = 0; i < lsize; ++i)
    {    
        ps_buffer[i] = ps[i];
        ps[i] = ps[i] + ps_dt[i]*dt;
        double pi = ps[i] - ptop; 

        // // Radiation
        // double albedo = 0.11;
        // albedo -= (terrain[i] > 0)*0.025;
        // albedo += (T[i + (nlevel - 2)*lsize] < 263)*0.2;
        // albedo = 1.4*((terrain[i] > 0) && (T[i + (nlevel - 2)*lsize] < 263))*albedo + 
        //         !(((terrain[i] > 0) && (T[i + (nlevel - 2)*lsize] < 263)))*albedo;

        // double SZA = sin(latitude[i])*sin(declination) + cos(latitude[i])*cos(declination)*cos(longitude[i] - hour_angle);
        // double sunlit = SZA > 0;
        // double rad_flux = dt*sunlit*(1 - albedo)*irradiance*SZA*g/(cp*ds);

        // // Evaporation
        double p = pi*(1 - ds/2) + ptop;
        // double q_sat = (pi/p)*610.94*exp(17.625*(T[i + (nlevel - 2)*lsize] - 273.15)/(T[i + (nlevel - 2)*lsize] - 30.08));
        // double wet = terrain[i] < 1 && q_sat > q[i + (nlevel - 2)*lsize];
        // double windspeed = sqrt(pow((u[i + (nlevel - 2)*lsize] + v[i + (nlevel - 2)*lsize])/pi, 2));
        // double evap_rate = (25 + 19*windspeed)/(60*60);
        // double evap = wet*evap_rate*dt*g*(q_sat - q[i + (nlevel - 2)*lsize])/(pi*ds);

        for (int level = 0; level < nlevel - 1; ++level)
        {
            int l_ = level*lsize;
            u_buffer[i + l_] = u[i + l_];
            v_buffer[i + l_] = v[i + l_];
            u[i + l_] = u[i + l_] + u_accel[i + l_]*dt;
            v[i + l_] = v[i + l_] + v_accel[i + l_]*dt;

            p = ps[i]*(sigma[level] + ds/2) + ptop;
            theta_buffer[i + l_] = theta[i + l_];
            // q_buffer[i + l_] = q[i + l_];
            theta[i + l_] = theta[i + l_] + theta_dt[i + l_]*dt;
            // q[i + l_] = q[i + l_] + q_dt[i + l_]*dt;
            T[i + l_] = (theta[i + l_])*pow(p_ref/p, -R/cp);

            // // Condensation
            // p = ps[i]*(sigma[level] + ds/2) + ptop;
            // q_sat = (pi/p)*610.94*exp(17.625*(T[i + l_] - 273.15)/(T[i + l_] - 30.08));
            // double condens_rate = q[i + l_] > q_sat ? 10 : 0;
            // double condens = dt*condens_rate*(q[i + l_] - q_sat);
            // double ground = level == nlevel - 2;
            // evap = 0;
            // condens = 0;
            // q[i + l_] += ground*evap - condens;

            // double T_ = T[i + l_] + (ground*rad_flux + LH*(condens - ground*evap))/pi;
            // T[i + l_] = T[i + l_]/pow(1 + 3*boltz*dt*GH*pow(T_, 3)*g/(cp*ds), 1/3);
            // theta[i + l_] = pi*T[i + l_]*pow(p_ref/p, R/cp);
        }
    }
}

void SkyCubeMP::LeapFrog(double dt)
{
    int tsize = nlat*nlat;
    int lsize = 6*tsize;

    #pragma omp for
    for (int i = 0; i < lsize; ++i)
    {   
        int pos = i%tsize;
        double p_swap = ps[i];
        double pi = ps[i] - ptop;
        ps[i] = ps_buffer[i] + ps_dt[i]*(2*dt); 
        ps_buffer[i] = p_swap + timefilter*(ps[i] -2*p_swap + ps_buffer[i]);

        // // Radiation
        // double albedo = 0.11;
        // albedo -= (terrain[i] > 0)*0.025;
        // albedo += (T[i + (nlevel - 2)*lsize] < 263)*0.2;
        // albedo = 1.4*((terrain[i] > 0) && (T[i + (nlevel - 2)*lsize] < 263))*albedo + 
        //         !(((terrain[i] > 0) && (T[i + (nlevel - 2)*lsize] < 263)))*albedo;

        // double SZA = sin(latitude[i])*sin(declination) + cos(latitude[i])*cos(declination)*cos(longitude[i] - hour_angle);
        // double sunlit = SZA > 0;
        // double rad_flux = dt*sunlit*(1 - albedo)*irradiance*SZA*g/(cp*ds);

        // // Evaporation
        double p = ps[i]*(1 - ds/2) + ptop;
        // double q_sat = (pi/p)*610.94*exp(17.625*(T[i + (nlevel - 2)*lsize] - 273.15)/(T[i + (nlevel - 2)*lsize] - 30.08));
        // double wet = terrain[i] < 1 && q_sat > q[i + (nlevel - 2)*lsize];
        // double windspeed = sqrt(pow((u[i + (nlevel - 2)*lsize] + v[i + (nlevel - 2)*lsize])/pi, 2));
        // double evap_rate = (wet*g/(pi*ds))*(25 + 19*windspeed)/(60*60);

        for (int level = 0; level < nlevel - 1; ++level)
        {
            int l_ = level*lsize;
            double u_swap = u[i + l_];
            double v_swap = v[i + l_];
            double theta_swap = theta[i + l_];
            double q_swap = q[i + l_];
            u[i + l_] = u_buffer[i + l_] + u_accel[i + l_]*(2*dt);
            v[i + l_] = v_buffer[i + l_] + v_accel[i + l_]*(2*dt);

            double pi = ps[i] - ptop;
            double p = ps[i]*(sigma[level] + ds/2) + ptop;
            theta[i + l_] = theta_buffer[i + l_] + theta_dt[i + l_]*(2*dt);

            // implicit schemes
            int north_halo = 1 - 2*(pos%nlat == 0);
            int south_halo = 1 - 2*(pos%nlat == nlat - 1);
            int east_halo = 1 - 2*(pos/nlat == nlat - 1);
            int west_halo = 1 - 2*(pos/nlat == 0);
            int northeast_halo = (pos != nlat*(nlat - 1))*north_halo;
            int southwest_halo = (pos != nlat - 1)*south_halo;
            int plate = (i%lsize)/tsize;
            int axis_flip = 1 - 2*(plate == 2 || plate == 3 || plate == 5);
            
            double f = 2*OM*axis_flip*sin(latitude[i]);

            double u_interp = (u[i + l_] + west_halo*(*(uwest[i] + l_))
            + south_halo*(*(usouth[i] + l_)) + southwest_halo*(*(usouthwest[i] + l_)))/4;

            double v_interp = (v[i + l_] + north_halo*(*(vnorth[i] + l_))
            + east_halo*(*(veast[i] + l_)) + northeast_halo*(*(vnortheast[i] + l_)))/4;

            u[i + l_] = (u[i + l_] + f*dt*v_interp)/(pow(f*dt, 2) + friction[i + l_]*dt + 1);
            v[i + l_] = (v[i + l_] - f*dt*u_interp)/(pow(f*dt, 2) + friction[i + l_]*dt + 1);

            // // Condensation
            // p = ps[i]*(sigma[level] + ds/2) + ptop;
            // q_sat = (pi/p)*610.94*exp(17.625*(T[i + l_] - 273.15)/(T[i + l_] - 30.08));
            // double condens_rate = q[i + l_] > q_sat ? 10 : 0;
            // double ground = level == nlevel - 2;
            // double& old_q = q[i + l_];
            // q[i + l_] = (q[i + l_] + dt*(ground*evap_rate - condens_rate)*q_sat)
            //                          /(1 + dt*(ground*evap_rate + condens_rate));

            // double T_ = T[i + l_] + (ground*rad_flux + LH*(old_q - q[i + l_]))/pi;
            // T[i + l_] = T_/cbrt(1 + 3*boltz*dt*GH*pow(T_, 3)*g/(cp*pi*ds));
            // theta[i + l_] = pi*T[i + l_]*pow(p_ref/p, R/cp);

            // buffer swap
            u_buffer[i + l_] = u_swap + timefilter*(u[i + l_] -2*u_swap + u_buffer[i + l_]);
            v_buffer[i + l_] = v_swap + timefilter*(v[i + l_] -2*v_swap + v_buffer[i + l_]);
            theta_buffer[i + l_] = theta_swap + timefilter*(theta[i + l_] -2*theta_swap + theta_buffer[i + l_]);
            q_buffer[i + l_] = q_swap + timefilter*(q[i + l_] -2*q_swap + q_buffer[i + l_]);
        }
    }
}

void SkyCubeMP::Print(int face)
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
            }
        }
    }
}

void SkyCubeMP::WriteState()
{
    std::ofstream ufile("u.csv", std::ios_base::app);
    std::ofstream vfile("v.csv", std::ios_base::app);
    std::ofstream sfile("s.csv", std::ios_base::app);
    std::ofstream pfile("p.csv", std::ios_base::app);
    std::ofstream gfile("g.csv", std::ios_base::app);
    std::ofstream tfile("t.csv", std::ios_base::app);
    std::ofstream qfile("q.csv", std::ios_base::app);
    std::ostringstream streamObj;
    std::string strObj;
    
    int tsize = nlat*nlat;
    int lsize = 6*tsize;

    for (int level = 0; level < nlevel - 1; ++level)
    {
        int l_ = level*lsize;

        if (level < nlevel/2 - 1)
        {
            continue;
        }

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

            streamObj << std::scientific << std::setprecision(8);
            streamObj << geopot[i + l_];
            // streamObj << (radius/100000)*terrain[i];
            strObj = streamObj.str();
            gfile << strObj << ',';
            streamObj.str(std::string());

            streamObj << std::scientific << std::setprecision(5);
            streamObj << theta[i + l_];
            strObj = streamObj.str();
            tfile << strObj << ',';
            streamObj.str(std::string());

            streamObj << std::scientific << std::setprecision(5);
            streamObj << q[i + l_]/pi;
            strObj = streamObj.str();
            qfile << strObj << ',';
            streamObj.str(std::string());
        }

        ufile << std::endl;
        vfile << std::endl;
        sfile << std::endl;
        gfile << std::endl;
        tfile << std::endl;
        qfile << std::endl;
    }

    streamObj << std::scientific << std::setprecision(8);

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
    qfile << std::endl;

    ufile.close();
    vfile.close();
    sfile.close();
    gfile.close();
    pfile.close();
    tfile.close();
    qfile.close();
}

SkyCubeMP::SkyCubeMP(int nl, int n_lev, double rad, double p0, double pt, std::vector<int>& terre)
:
nlat(nl),
nlevel(n_lev),
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
q(new double[6*nlat*nlat*(nlevel - 1)]),
q_dt(new double[6*nlat*nlat*(nlevel - 1)]),
q_buffer(new double[6*nlat*nlat*(nlevel - 1)]),
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

    double Ts = 298;

    for (int i = 0; i < lsize; ++i)
    {
        int l_ = lsize*(nlevel - 1);

        sigma_vel[i + l_] = 0;
        surface_pot[i] = (radius/10)*terrain[i]*g;
        ps[i] = pow(pow(p_ref, R/cp)*(1 - surface_pot[i]/(cp*Ts)), cp/R);
        ps_buffer[i] = ps[i];
        ps_dt[i] = 0;
        q[i + l_] = 0;
        q_buffer[i + l_] = 0;
        q_dt[i + l_] = 0;
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
                Wx2j[pos] = &y1i[iwest[pos]%tsize];
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
                Ny1j[pos] = &x2i[inorth[pos]%tsize];

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

    for (int i = 0; i < lsize; ++i)
    {
        double pi = ps[i] - ptop;

        for (int level = 0; level < nlevel - 1; ++level)
        {
            int l_ = 6*nlat*nlat*level;

            theta[i + l_] = Ts;
            theta_buffer[i + l_] = theta[i + l_];

            double p = (sigma[level] + ds/2)*pi + ptop;
            T[i + l_] = theta[i + l_]*pow((p_ref/p), -R/cp);
            friction[i + l_] = level < nlevel - 2 ? 0 : (1 + terrain[i]/10)*0.000015;
        }
    }

    // ps[nlat*(nlat/2) + nlat/2] = ps[nlat*(nlat/2) + nlat/2] + 100;
    // ps_buffer[nlat*(nlat/2) + nlat/2] = ps_buffer[nlat*(nlat/2) + nlat/2] + 100;
    theta[nlat*(nlat/2) + nlat/2] = theta[nlat*(nlat/2) + nlat/2] + 20;
    theta_buffer[nlat*(nlat/2) + nlat/2] = theta_buffer[nlat*(nlat/2) + nlat/2] + 20;
}

SkyCubeMP::~SkyCubeMP()
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
    delete q;
    delete q_dt;
    delete q_buffer;
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