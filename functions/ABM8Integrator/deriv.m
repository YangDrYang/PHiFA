%------------------------------------------------------------------------------
%
% deriv.m :Computes the derivative of the state vector for the normalized
%          (GM=1) Kepler's problem in three dimensions
%
% Inputs:
%   t : time(s)
%   y : state vector(x,y,z,vx,vy,vz)
%
% Outputs:
%   yp : derivative of the state vector(vx,vy,vz,ax,ay,az)
%
% Last modified:   2018/01/27   M. Mahooti
%------------------------------------------------------------------------------
function yp = deriv(t, y)

fprintf('### Accel: time %f\n', t);

% State vector derivative
r = y(1:3);
v = y(4:6);
yp = [v;-r/((norm(r))^3)];

