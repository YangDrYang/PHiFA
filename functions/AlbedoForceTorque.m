%
%  albedo_calc.c
%  D-SPOSE
%
%  Created by Luc Sagnieres on 2018-07-02.
%  Copyright © 2018 Luc Sagnieres. All rights reserved.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FUNCTION NAME:        albedo_calc.c
%
% DESCRIPTION:          This function calculates the acceleration and
%                       torque due to reflected and emitted Earth radiation
%                       (See Section 2.3.3 in Sagnieres (2018) Doctoral Thesis)
%
% AUTHOR:               Luc Sagnieres
% DATE:                 July 2, 2018
% VERSION:              1
% AUTHOR:               Yang Yang
% DATE:                 July 10, 2019
% VERSION:              2 C->Matlab
%
% INPUT:                double p[3]: 3x1 position vector (m)
%                       double mjd[1]: day
%                       double r_sun[3]: position of Sun (m)
%                       double C_ecef2teme[3][3]: rotation matrix from ECEF
%                         frame to TEME
%                       double C_i2b[3][3]: rotation matrix from inertial
%                         frame to body frame
%                       int n_surf: number of surfaces in geometry model
%                       struct surface geometry[n_surf]: surface geometry model
%                       double m: mass (kg)
%                       double albedo[12][20][40][2]: albedo and IR coefficients
%                       int in_alb_a: inclusion of albedo acceleration?
%                       int in_alb_g: inclusion of albedo torque?
%                       int in_ir_a: inclusion of IR acceleration?
%                       int in_ir_g: inclusion of IR torque?
%
% OUTPUT:               double a_alb[3]: reflected Earth radiation acceleration
%                         in inertial frame
%                       double g_alb[3]: reflected Earth radiation torque
%                         in body-fixed frame
%                       double a_in[3]: emitted Earth radiation acceleration
%                         in inertial frame
%                       double g_in[3]: emitted Earth radiation torque
%                         in body-fixed frame
%                       double proA_alb[1]: projected area for reflected Earth 
%                          radiation
%                       double proA_in[1]: projected area for emitted Earth 
%                          radiation
%
% COUPLING:             - surface.h
%                       - vectors2angle.c
%                       - norm.c
%                       - crossprod.c
%                       - matxvec.c
%                       - invertmat.c
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [a_alb,g_alb,a_ir,g_ir,proA_alb,proA_ir] = AlbedoForceTorque(obj,p, mjd, r_sun, ...
C_ecef2teme, C_i2b, n_surf, facet, m, albedo, in_alb_a, in_alb_g, in_ir_a, in_ir_g)

a_alb = zeros(3,1);
g_alb = zeros(3,1);
a_ir = zeros(3,1);
g_ir = zeros(3,1);
proA_alb=0;
proA_ir=0;
C_b2i = C_i2b';
area_ir = 0;
area_alb = 0;

% Initialize acceleration for each surface
f_surf = zeros(n_surf,2,3);

% Surface area for each Grid (9° x 9°) as a function of latitude (90-81 to 9-0) independent of longitude (in m^2)
grid_area = [79189845238,235467450694,385524392958,525488995462,651870341413,...
    761652997845,852348921372,922010908066,969218205729,993048212774,...
    993048212774,969218205729,922010908066,852348921372,761652997845,...
    651870341413,525488995462,385524392958,235467450694,79189845238];

% Earth Equatorial Radius and Eccentricity
re         =     6378137;
eesqrd     =     0.006694385000;

% Interpolate albedo to current day
[yy,mm,dd,hh,min,ss] = invjday(mjd+2400000.5);
doy = floor(finddays(yy,mm,dd,hh,min,ss));
albedo_interp = zeros(20,40,2);

for i=1:20
    for j=1:40
        if (doy<=15)
            albedo_interp(i,j,1) = albedo(12,i,j,1) + (doy + 15) * (albedo(1,i,j,1) - albedo(12,i,j,1)) / 31.0;
            albedo_interp(i,j,2) = albedo(12,i,j,2) + (doy + 15) * (albedo(1,i,j,2) - albedo(12,i,j,2)) / 31.0;
        elseif (doy<=46)
            albedo_interp(i,j,1) = albedo(1,i,j,1) + (doy - 15) * (albedo(2,i,j,1) - albedo(1,i,j,1)) / 31.0;
            albedo_interp(i,j,2) = albedo(1,i,j,2) + (doy - 15) * (albedo(2,i,j,2) - albedo(1,i,j,2)) / 31.0;
        elseif (doy<=75)
            albedo_interp(i,j,1) = albedo(2,i,j,1) + (doy - 46) * (albedo(3,i,j,1) - albedo(2,i,j,1)) / 29.0;
            albedo_interp(i,j,2) = albedo(2,i,j,2) + (doy - 46) * (albedo(3,i,j,1) - albedo(2,i,j,2)) / 29.0;
        elseif (doy<=106)
            albedo_interp(i,j,1) = albedo(3,i,j,1) + (doy - 75) * (albedo(4,i,j,1) - albedo(3,i,j,1)) / 31.0;
            albedo_interp(i,j,2) = albedo(3,i,j,2) + (doy - 75) * (albedo(4,i,j,2) - albedo(3,i,j,2)) / 31.0;
        elseif (doy<=136)
            albedo_interp(i,j,1) = albedo(4,i,j,1) + (doy - 106) * (albedo(5,i,j,1) - albedo(4,i,j,1)) / 30.0;
            albedo_interp(i,j,2) = albedo(4,i,j,2) + (doy - 106) * (albedo(5,i,j,2) - albedo(4,i,j,2)) / 30.0;
        elseif (doy<=167)
            albedo_interp(i,j,1) = albedo(5,i,j,1) + (doy - 136) * (albedo(6,i,j,1) - albedo(5,i,j,1)) / 31.0;
            albedo_interp(i,j,2) = albedo(5,i,j,2) + (doy - 136) * (albedo(6,i,j,2) - albedo(5,i,j,2)) / 31.0;
        elseif (doy<=197)
            albedo_interp(i,j,1) = albedo(6,i,j,1) + (doy - 167) * (albedo(7,i,j,1) - albedo(6,i,j,1)) / 30.0;
            albedo_interp(i,j,2) = albedo(6,i,j,2) + (doy - 167) * (albedo(7,i,j,2) - albedo(6,i,j,2)) / 30.0;            
        elseif (doy<=228)
            albedo_interp(i,j,1) = albedo(7,i,j,1) + (doy - 197) * (albedo(8,i,j,1) - albedo(7,i,j,1)) / 31.0;
            albedo_interp(i,j,2) = albedo(7,i,j,2) + (doy - 197) * (albedo(8,i,j,2) - albedo(7,i,j,2)) / 31.0;
        elseif (doy<=259)
            albedo_interp(i,j,1) = albedo(8,i,j,1) + (doy - 228) * (albedo(9,i,j,1) - albedo(8,i,j,1)) / 31.0;
            albedo_interp(i,j,2) = albedo(8,i,j,2) + (doy - 228) * (albedo(9,i,j,2) - albedo(8,i,j,2)) / 31.0;
        elseif (doy<=289)
            albedo_interp(i,j,1) = albedo(9,i,j,1) + (doy - 259) * (albedo(10,i,j,1) - albedo(9,i,j,1)) / 30.0;
            albedo_interp(i,j,2) = albedo(9,i,j,2) + (doy - 259) * (albedo(10,i,j,2) - albedo(9,i,j,2)) / 30.0;
        elseif (doy<=320)
            albedo_interp(i,j,1) = albedo(10,i,j,1) + (doy - 289) * (albedo(11,i,j,1) - albedo(10,i,j,1)) / 31.0;
            albedo_interp(i,j,2) = albedo(10,i,j,2) + (doy - 289) * (albedo(11,i,j,2) - albedo(10,i,j,2)) / 31.0;
        elseif (doy<=350)
            albedo_interp(i,j,1) = albedo(11,i,j,1) + (doy - 320) * (albedo(12,i,j,1) - albedo(11,i,j,1)) / 30.0;
            albedo_interp(i,j,2) = albedo(11,i,j,2) + (doy - 320) * (albedo(12,i,j,2) - albedo(11,i,j,2)) / 30.0;            
        else
            albedo_interp(i,j,1) = albedo(12,i,j,1) + (doy - 350) * (albedo(1,i,j,1) - albedo(12,i,j,1)) / 31.0;
            albedo_interp(i,j,2) = albedo(12,i,j,2) + (doy - 350) * (albedo(1,i,j,2) - albedo(12,i,j,2)) / 31.0;            
        end
    end
end

% Calculate grid center position in ECEF (approximate outward surface normal as same divided by radius)
grid_position_ecef = zeros(20,40,3);
grid_position_eci = zeros(20,40,3);
norm_grid_position = zeros(20,40);
unit_grid_eci = zeros(20,40,3);
for i=1:20
    for j=1:40
        grid_position_ecef(i,j,1) = re*cos((90-4.5-9*(i-1))*pi/180.0)*cos((4.5+9*(j-1))*pi/180.0);
        grid_position_ecef(i,j,2) = re*cos((90-4.5-9*(i-1))*pi/180.0)*sin((4.5+9*(j-1))*pi/180.0);
        grid_position_ecef(i,j,3) = re*sqrt(1-eesqrd)*sin((90-4.5-9*(i-1))*pi/180.0);
        
        tmp(1:3,1) = grid_position_ecef(i,j,:);
        norm_grid_position(i,j) = norm(tmp);
        grid_position_eci(i,j,:) = C_ecef2teme*tmp;

        unit_grid_eci(i,j,:) = grid_position_eci(i,j,:)/norm_grid_position(i,j);
    end
end

% Solar Flux
rsun = norm(r_sun);
phi = 1361*(149597870700.0/rsun)*(149597870700.0/rsun);

% For each Earth grid
for i=1:20
    for j=1:40

        % Calculate grid-satellite distance
        grid_sat_pos(1:3,1) = p - reshape(grid_position_eci(i,j,:),3,1);
        grid_sat_dist = norm(grid_sat_pos);

        % Check if grid is in view of satellite
        grid_view_sat = vectors2angle(reshape(unit_grid_eci(i,j,:),3,1),grid_sat_pos);
        cos_grid_view_sat = cos(grid_view_sat);

        % If grid is in view of satellite
        if (cos_grid_view_sat > 0)

            % Check if Earth grid is in view of sun
            grid_view_sun = vectors2angle(reshape(unit_grid_eci(i,j,:),3,1),r_sun);
            cos_grid_view_sun = cos(grid_view_sun);

            % For every satellite surface
            for k=1:n_surf

                % Rotate surface unit normal in ECI frame
                surf_unit_eci = C_b2i*facet(k).normal;

                % Check if surface is in view of Earth grid
                grid_view_surf = vectors2angle(reshape(unit_grid_eci(i,j,:),3,1),surf_unit_eci);
                cos_grid_view_surf = cos(grid_view_surf);

                % If surface is in view of Earth grid
                if (cos_grid_view_surf>0)

                    % Optical properties
                    crd = obj.surfAtr.alpha*(1-obj.surfAtr.beta);
                    crs = obj.surfAtr.alpha*obj.surfAtr.beta;
                    ca = (1-obj.surfAtr.alpha);
                    
                    % Optical properties
                    crd_ir = obj.surfAtr.alpha*(1-obj.surfAtr.beta);
                    crs_ir = obj.surfAtr.alpha*obj.surfAtr.beta;
                    ca_ir = (1-obj.surfAtr.alpha);

                    % Speed of light
                    c = 299792458.0;

                    % Projected surface area
                    Ap = facet.area*cos_grid_view_surf;

                    % Calculate IR acceleration for that surface from that surface grid in ECI
                    if (in_ir_a || in_ir_g)

                        % Calculate flux
                        flux_ir = phi*albedo_interp(i,j,2)/1000.0*cos_grid_view_sat*grid_area(i)/(4*pi*grid_sat_dist*grid_sat_dist);

                        % Calculate force on that surface
                        for ii = 1:3
                            f_surf(k,2,ii) = f_surf(k,2,ii) + flux_ir/c*Ap*((ca_ir+crd_ir)*unit_grid_eci(i,j,ii)...
                                + (2*crd_ir/3.0 + 2*crs_ir*cos_grid_view_surf)*surf_unit_eci(ii));
                        end
                        if area_ir < Ap
                            area_ir = Ap;
                        end
                    end

                    % If Earth grid is in view of sun
                    if (cos_grid_view_sun>0)

                        % Calculate albedo acceleration for that surface from that surface grid in ECI
                        if (in_alb_a || in_alb_g)

                            % Calculate flux
                            flux_alb = phi*albedo_interp(i,j,1)/1000.0*cos_grid_view_sat*cos_grid_view_sun*grid_area(i)/(pi*grid_sat_dist*grid_sat_dist);

                            % Calculate force on that surface
                            for ii = 1:3
                                f_surf(k,1,ii) = f_surf(k,1,ii) + flux_alb/c *Ap*((ca+crd)*unit_grid_eci(i,j,ii)...
                                + (2*crd/3.0 + 2*crs*cos_grid_view_surf)*surf_unit_eci(ii));
                            end
                            
                            if area_alb < Ap
                                area_alb = Ap;
                            end
                        end
                    end
                end
            end
        end
    end
end

% Calculate total IR acceleration
if (in_ir_a)
    for k=1:n_surf
        a_ir = a_ir+reshape(f_surf(k,2,:),3,1)/m;
    end
    proA_ir = area_ir;
end

% Calculate total albedo acceleration
if (in_alb_a)
    for k=1:n_surf
        a_alb = a_alb+reshape(f_surf(k,1,:),3,1)/m;
    end
    proA_alb = area_alb;
end

% % Calculate center of pressure
% cp_b = zeros(n_surf,3);
% for i=1:n_surf
%     cp_b(i,:) = sum(geometry(i).vertices(0,:))/3.0;
% end

% Calculate IR torque for each surface in BODY
if (in_ir_g)
    for k=1:n_surf
        % Force in body frame
        f_surf_body = C_i2b*reshape(f_surf(k,2,:),3,1);

        % Calculate torque
        g_surf = cross(facet.areabarycenter+obj.offset,f_surf_body);
        g_ir = g_ir+g_surf;
    end
end

% Calculate albedo torque for each surface in BODY
if (in_alb_g)
    for k=1:n_surf
        % Force in body frame
        f_surf_body = C_i2b*reshape(f_surf(k,1,:),3,1);

        % Calculate torque
        g_surf = cross(facet.areabarycenter+obj.offset,f_surf_body);
        g_alb = g_alb+g_surf;
    end
end
