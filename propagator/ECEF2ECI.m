%--------------------------------------------------------------------------
%
% ECEF2ECI: Transforms Earth Centered Earth Fixed (ECEF) coordinates to 
%           Earth Centered Inertial (ECI) coordinates
%
% Last modified:   2018/01/27   M. Mahooti
%
%--------------------------------------------------------------------------
function Y = ECEF2ECI(MJD_UTC, Y0)
[U,dU] = transMatEci2Ecef(MJD_UTC);

% Transformation from WGS to ICRS
r = U'*Y0(1:3)';
v = U'*Y0(4:6)' + dU'*Y0(1:3)';
Y = [r;v];

