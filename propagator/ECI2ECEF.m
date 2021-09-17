%--------------------------------------------------------------------------
%
% ECI2ECEF: Transforms Earth Centered Inertial (ECI) coordinates to Earth
%           Centered Earth Fixed (ECEF) coordinates
%
% Last modified:   2018/01/27   M. Mahooti
%
%--------------------------------------------------------------------------
function Y = ECI2ECEF(MJD_UTC, Y0)
[U,dU] = transMatEci2Ecef(MJD_UTC);

% Transformation from ICRS to WGS
r = U*Y0(1:3)';
v = U*Y0(4:6)' + dU*Y0(1:3)';
Y = [r;v];

