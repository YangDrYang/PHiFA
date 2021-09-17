% Constants
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FUNCTION NAME:        magnetic_field.c
%
% DESCRIPTION:          This function will calculate the magnetic field
%                       vector and its time derivative at the satellite position
%
% AUTHOR:               Luc Sagnieres
% DATE:                 September 11, 2016
% VERSION:              1
% AUTHOR:               Yang Yang
% DARE:                 July 10, 2019
% VWRSION:              2 C->Matlab
%
% INPUT:                double mJD: modified Julian date
%                       double p_LLA[3]: position in ECEF frame (m)
%                       double v_LLA[3]: velocity in ECEF frame (m s-1)
%                       double LLA[4]:
%                         - LLA[0]: geocentric latitude (rad)
%                         - LLA[1]: geodetic latitude (rad)
%                         - LLA[2]: longitude (rad)
%                         - LLA[3]: altitude (km)
%                       double C_ecef2teme[3][3]: rotation matrix from ECEF to TEME
%                       double G[14][14][25]: magnetic potential coefficients
%                       double H[14][14][25]: magnetic potential coefficients
%
% OUTPUT:               double B_field_i[3]: magnetic field vector in inertial frame
%                       double B_field_i_dot[3]: time derivative of magnetic field
%                         vector in inertial frame as seen from orbiting
%                         spacecraft


function [B_field_i,B_field_dot_i] = magnet_field(mJD, p_LLA, v_LLA, LLA, C_ecef2teme, G, H)

a = 6371200; % Geomagnetic conventional Earth's mean reference spherical radius

% Get Decimal Year

year_decimal = decyear(JD2Date(mJD+2400000.5));

% Interpolate G and H linearly
G_t = zeros(14,14);
H_t = zeros(14,14);
if (floor(year_decimal) < 2005)
    k = 20;
elseif (floor(year_decimal) < 2010)
    k = 21;
elseif (floor(year_decimal) < 2015)
    k = 22;
else
    k = 23;
end
year_start = (k+380)*5;

if (k<23)
    for i=1:14
        for j=1:14
            G_t(i,j) = G(i,j,k+1)+(year_decimal-year_start)*(G(i,j,k+2)-G(i,j,k+1))/5.0;
            H_t(i,j) = H(i,j,k+1)+(year_decimal-year_start)*(H(i,j,k+2)-H(i,j,k+1))/5.0;
        end
    end
elseif (k==23)
    for i=1:14
        for j=1:14
            G_t(i,j) = G(i,j,k+1)+(year_decimal-year_start)*G(i,j,k+2);
            H_t(i,j) = H(i,j,k+1)+(year_decimal-year_start)*H(i,j,k+2);

        end
    end
end

%Time derivative of G and H
dG_t = zeros(14,14);
dH_t = zeros(14,14);
if (k<23)
    for i=1:14
        for j=1:14
            dG_t(i,j) = (G(i,j,k+2)-G(i,j,k+1))/(5.0*365.25*24*60*60);
            dH_t(i,j) = (H(i,j,k+2)-H(i,j,k+1))/(5.0*365.25*24*60*60);
        end
    end
elseif (k==23)
    for i=1:14
        for j=1:14
            dG_t(i,j) = G(i,j,k+2)/(365.25*24*60*60);
            dH_t(i,j) = H(i,j,k+2)/(365.25*24*60*60);
        end
    end
end

%Spherical geocentric distance, longitude and latitude (declination), and colatitude (m and rad)
r = norm(p_LLA);
lon = LLA(3);
lat = asin(p_LLA(3)/r);
colat = pi/2.0 - lat;

% Check for Poles
if (colat>-0.00000001*pi/180.0 && colat<0.00000001*pi/180.0)
    colat=0.00000001*pi/180.0;
elseif(colat<180.00000001*pi/180.0 && colat>179.99999999*pi/180.0)
    colat=179.99999999*pi/180.0;
end

%Calculate Schmidt, associated Legendre, quasi-normalized blablabla
P = zeros(14,14);
P(1,1) = 1;
P(2,1) = cos(colat);
for n=2:14
    P(n,n) = sin(colat)*P(n-1,n-1);
end
for n=3:14
    for m=1:n-1
        P(n,m) = cos(colat)*P(n-1,m) - ((n-2)*(n-2)-(m-1)*(m-1))*P(n-2,m)/((2*(n-1)-1)*(2*(n-1)-3));
    end
end

dP = zeros(14,14);
dP(1,1) = 0;
dP(2,1) = -sin(colat);
for n=2:14
    dP(n,n) = sin(colat)*dP(n-1,n-1) + cos(colat)*P(n-1,n-1);
end
for n=3:14
    for m=1:n-1
        dP(n,m) = cos(colat)*dP(n-1,m) - sin(colat)*P(n-1,m) - ((n-2)*(n-2)-(m-1)*(m-1))*dP(n-2,m)/((2*(n-1)-1)*(2*(n-1)-3));
    end
end

ddP = zeros(14,14);
ddP(1,1) = 0;
ddP(2,1) = -cos(colat);
for n=2:14
    ddP(n,n) = 2*cos(colat)*dP(n-1,n-1) - sin(colat)*P(n-1,n-1) + sin(colat)*ddP(n-1,n-1);
end
for n=3:14
    for m=1:n-1
        ddP(n,m) = cos(colat)*ddP(n-1,m) - 2*sin(colat)*dP(n-1,m) - cos(colat)*P(n-1,m) - ((n-2)*(n-2)-(m-1)*(m-1))*ddP(n-2,m)/((2*(n-1)-1)*(2*(n-1)-3));
    end
end
S = zeros(14,14);
S(1,1) = 1;
for n=2:14
    S(n,1) = S(n-1,1)*(2*(n-1)-1)/(n-1);
end
for n=2:14
    for m=2:n
        if (m==2)
            delta = 1;
        else
            delta = 0;
        end
        S(n,m) = S(n,m-1)*sqrt(((n-1)-(m-1)+1)*(delta+1)/((n-1)+(m-1)));
    end
end

for n=1:14
    for m=1:n
        P(n,m) = P(n,m) *S(n,m);
        dP(n,m) = dP(n,m) *S(n,m);
        ddP(n,m) = ddP(n,m) *S(n,m);
    end
end

% Partial Derivatives
dVdr = 0;
dVd0 = 0;
dVdl = 0;
for l = 2:14
    for m=1:l
        dVdr = dVdr - (a/r)^(l+1)*l*P(l,m)*(G_t(l,m)*cos((m-1)*lon)+H_t(l,m)*sin((m-1)*lon));
        dVd0 = dVd0 + a*(a/r)^l*dP(l,m)*(G_t(l,m)*cos((m-1)*lon)+H_t(l,m)*sin((m-1)*lon));
        dVdl = dVdl + a*(a/r)^l*(m-1)*P(l,m)*(H_t(l,m)*cos((m-1)*lon)-G_t(l,m)*sin((m-1)*lon));
    end
end

% Second Partial Derivatives
d2Vdr2 = 0;
d2Vd0dr = 0;
d2Vdldr = 0;
d2Vdtdr = 0;
d2Vd02 = 0;
d2Vdld0 = 0;
d2Vdtd0 = 0;
d2Vdl2 = 0;
d2Vdtdl = 0;
for l = 2:14
    for m = 1:l
        d2Vdr2 = d2Vdr2 + 1/r*(a/r)^(l+1)*(l+1)*l*P(l,m)*(G_t(l,m)*cos((m-1)*lon)+H_t(l,m)*sin((m-1)*lon));
        d2Vd0dr = d2Vd0dr - (a/r)^(l+1)*l*dP(l,m)*(G_t(l,m)*cos((m-1)*lon)+H_t(l,m)*sin((m-1)*lon));
        d2Vdldr = d2Vdldr + (a/r)^(l+1)*l*(m-1)*P(l,m)*(G_t(l,m)*sin((m-1)*lon)-H_t(l,m)*cos((m-1)*lon));
        d2Vdtdr = d2Vdtdr - (a/r)^(l+1)*l*P(l,m)*(dG_t(l,m)*cos((m-1)*lon)+dH_t(l,m)*sin((m-1)*lon));
        d2Vd02 = d2Vd02 + a*(a/r)^l*ddP(l,m)*(G_t(l,m)*cos((m-1)*lon)+H_t(l,m)*sin((m-1)*lon));
        d2Vdld0 = d2Vdld0 + a*(a/r)^l*(m-1)*dP(l,m)*(H_t(l,m)*cos((m-1)*lon)-G_t(l,m)*sin((m-1)*lon));
        d2Vdtd0 = d2Vdtd0 + a*(a/r)^l*dP(l,m)*(dG_t(l,m)*cos((m-1)*lon)+dH_t(l,m)*sin((m-1)*lon));
        d2Vdl2 = d2Vdl2 - a*(a/r)^l*(m-1)*(m-1)*P(l,m)*(G_t(l,m)*cos((m-1)*lon)+H_t(l,m)*sin((m-1)*lon));
        d2Vdtdl = d2Vdtdl + a*(a/r)^l*(m-1)*P(l,m)*(dH_t(l,m)*cos((m-1)*lon)-dG_t(l,m)*sin((m-1)*lon));
    end
end

%Magnetic Field in Spherical Components
B_field_spherical = zeros(3,1);
B_field_spherical(1) = -dVdr*1.0e-9;
B_field_spherical(2) = -dVd0/r*1.0e-9;
B_field_spherical(3) = -dVdl/(r*sin(colat))*1.0e-9;

%Magnetic Field in North, East, Down Components
B_field_ned = zeros(3,1);
B_field_ned(1) = -B_field_spherical(2);
B_field_ned(2) = B_field_spherical(3);
B_field_ned(3) = -B_field_spherical(1);

%Rotate to ecef frame
C_ned2ecef = ...
    [-sin(lat)*cos(lon),    -sin(lon),  -cos(lat)*cos(lon);...
    -sin(lat)*sin(lon),    cos(lon),   -cos(lat)*sin(lon);...
    cos(lat),              0,          -sin(lat)];
B_field_ecef = C_ned2ecef*B_field_ned;
B_field_i = C_ecef2teme*B_field_ecef;

%Velocity in Spherical Coordinates
v_sph = zeros(3,1);
v_ned = C_ned2ecef'*v_LLA;
v_sph(1) = - v_ned(3);
v_sph(2) = - v_ned(1);
v_sph(3) = v_ned(2);

%Time Derivative of Magnetic Field
B_field_sph_dot = zeros(3,1);
B_field_sph_dot(1) = -(v_sph(1)*d2Vdr2 + v_sph(2)/r*d2Vd0dr + ...
    v_sph(3)/(r*sin(colat))*d2Vdldr + d2Vdtdr)*1.0e-9;
B_field_sph_dot(2) = -(v_sph(1)/r*d2Vd0dr - v_sph(1)/(r*r)*dVd0 + ...
    v_sph(2)/(r*r)*d2Vd02 + v_sph(3)/(r*r*sin(colat))*d2Vdld0 + d2Vdtd0/r)*1.0e-9;
B_field_sph_dot(3) = -(v_sph(1)/(r*sin(colat))*d2Vdldr - v_sph(1)/(r*r*sin(colat))*dVdl - ...
    (v_sph(2)*cos(colat))/(r*r*sin(colat)*sin(colat))*dVdl + v_sph(2)/(r*r*sin(colat)*sin(colat))*d2Vdld0 + ...
    v_sph(3)/(r*r*sin(colat)*sin(colat))*d2Vdl2 + d2Vdtdl/(r*sin(colat)))*1.0e-9;

%Time Derivative Magnetic Field in North, East, Down Components
B_field_ned_dot = zeros(3,1);
B_field_ned_dot(1) = -B_field_sph_dot(2);
B_field_ned_dot(2) = B_field_sph_dot(3);
B_field_ned_dot(3) = -B_field_sph_dot(1);

%Rotate to ecef frame
w_ned = [v_ned(2)/r; -v_ned(1)/r; 0];
wxB_ned = cross(w_ned,B_field_ned);
B_field_ned_dot = B_field_ned_dot + wxB_ned;
B_field_ecef_dot = C_ned2ecef*B_field_ned_dot;

%Rotate to inertial frame
omegaearth = [0; 0; 7.29211514670698e-05];
wxB_ecef = cross(omegaearth,B_field_ecef);
B_field_ecef_dot = B_field_ecef_dot + wxB_ecef;
B_field_dot_i = C_ecef2teme*B_field_ecef_dot;    