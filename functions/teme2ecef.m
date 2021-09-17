
%   teme2ecef.c
%   D-SPOSE
% 
%   Created by Luc Sagnieres on 2017-10-27.
%   Copyright © 2018 Luc Sagnieres. All rights reserved.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FUNCTION NAME:        teme2ecef.c
%
% DESCRIPTION:          This function converts the position and velocity
%                       vectors from the TEME frame to the ECEF frame
%                       taking into account the effects of sidereal time
%                       (GMST) and polar motion following Vallado (2013)
%                       Section 3.7
%
% AUTHOR:               Luc Sagnieres (following Vallado algorithm)
% DATE:                 October 27, 2017
% VERSION:              1
% VERSION:              2, c->matlab by Yang
%
% INPUT:                double r_teme[3]: position vector in TEME (m)
%                       double v_teme[3]: velocity vector in TEME (m s-1)
%                       double ttt: julian centuries of terrestrial time
%                       double t2000ut1: seconds since January 1, 2000, 00:00:00 UT1
%                       double xp: polar motion parameter from EOP
%                       double yp: polar motion parameter from EOP
%                       double lod: length of day from EOP
%
% OUTPUT:               double r_ecef[3]: position vector in ECEF (m)
%                       double v_ecef[3]: velocity vector in ECEF (m s-1)
%                       double C_ecef2teme[3][3]: rotation matrix from
%                         ECEF to TEME
%
% COUPLING:             - matxvec.c
%                       - transpose.c
%                       - crossprod.c
%                       - polarm.c
%                       - matrixmult.c
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [r_ecef,v_ecef,C_ecef2teme] = teme2ecef(r_teme, v_teme, ttt, t2000ut1, xp, yp, lod)
    
%      Calculate GMST from J2000+dUT1
    tut1 = (t2000ut1-(12*60*60.0))/(60*60*24*36525.0);
    temp = - 6.2e-6 * tut1 * tut1 * tut1 + 0.093104 * tut1 * tut1 + (876600.0 * 3600.0 + 8640184.812866) * tut1 + 67310.54841;
    temp = fmod( temp*pi/180.0/240.0, 2*pi );
    if (temp<0.0)
        temp = temp+2*pi;
    end
    gmst = temp;
    
%      Rotation matrix from PEF to TEME considering sidereal time with GMST
    st =  [cos(gmst), -sin(gmst), 0; sin(gmst), cos(gmst) 0; 0, 0, 1];
    st_t = st';
    
%      Rotation matrix from ECEF to PEF considering polar motion
    pm = polarm_(xp,yp,ttt,'80');
    pm_t = pm';
    
%      Earth angular velocity
    thetasa    = 7.29211514670698e-05 * (1.0  - lod/86400.0 );
    omegaearth = [0, 0, thetasa]';
    
%      Position transformation
    r_pef = st_t*r_teme;
    r_ecef = pm_t*r_pef;
    C_ecef2teme = st*pm;
    
%      Velocity transformation
    v_pef = st_t*v_teme;
    v_temp = cross(omegaearth,r_pef);
    v_pef = v_pef - v_temp;
    v_ecef = pm_t*v_pef;
    
end