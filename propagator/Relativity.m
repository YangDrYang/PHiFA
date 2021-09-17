%--------------------------------------------------------------------------
%
% Relativisty: Computes the perturbational acceleration due to relativistic
%              effects
%
% Inputs:
%   r           Satellite position vector
%   v           Satellite velocity vector
% 
% Output:
%   a    		Acceleration (a=d^2r/dt^2)
%
% Last modified:   2018/01/27   M. Mahooti
%
%--------------------------------------------------------------------------
function a = Relativity(r, v)

p = clPropagator.instance();
% Relative position vector of satellite w.r.t. point mass 
r_Sat = norm(r);
v_Sat = norm(v);

% Acceleration 
a = p.const.GM_Earth/(p.const.c_light^2*r_Sat^3)*((4*p.const.GM_Earth/r_Sat-v_Sat^2)*r+4*dot(r,v)*v);

