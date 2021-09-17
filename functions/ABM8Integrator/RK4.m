%--------------------------------------------------------------------------
%
% Runge Kutta 4th order implementation
%
% Inputs:
%        func: function handle
%        t_0 : initial time
%        y_0 : initial state vector (x,y,z,vx,vy,vz)
%        h   : step size of integration
% 
% Output:
%        y   : state vector at t_0+h
% 
% Last modified:   2018/01/27   M. Mahooti
%--------------------------------------------------------------------------
function y = RK4(func, t_0, y_0, h)

k_1 = func(t_0    , y_0          );
k_2 = func(t_0+h/2, y_0+(h/2)*k_1);
k_3 = func(t_0+h/2, y_0+(h/2)*k_2);
k_4 = func(t_0+h  , y_0+    h*k_3);

y = y_0 + (h/6)*(k_1 + 2*k_2 + 2*k_3 + k_4);