%--------------------------------------------------------------------------
% 
% Adams-Bashforth-Moulton 8th order
% 
% Inputs:
%        func   : function handle
%        f_hist : a matrix which includes 8 rows and 7 columns including
%                 time(s) and state vector(x,y,z,vx,vy,vz)
%        h      : step size of integration
% 
% Output:
%        Y      : state vector(x,y,z,vx,vy,vz) at t_0+h
% 
% Last modified:   2018/01/27   M. Mahooti
%--------------------------------------------------------------------------
% function Y = ABM8(func, f_hist, h)
function Y = ABM8(func, h, f_hist)

bash = [434241.0, -1152169.0, 2183877.0, -2664477.0, 2102243.0, -1041723.0, 295767.0, -36799.0];
moul = [36799.0, 139849.0, -121797.0, 123133.0, -88547.0, 41499.0, -11351.0, 1375.0];
divisor = 1.0/120960.0;

% Calculate the predictor using the Adams-Bashforth formula 
Y = f_hist(8,2:end)' + h*divisor* ...
    ( bash(1)*func(f_hist(8,1),f_hist(8,2:end)') + bash(2)*func(f_hist(7,1),f_hist(7,2:end)') + ...
      bash(3)*func(f_hist(6,1),f_hist(6,2:end)') + bash(4)*func(f_hist(5,1),f_hist(5,2:end)') + ...
      bash(5)*func(f_hist(4,1),f_hist(4,2:end)') + bash(6)*func(f_hist(3,1),f_hist(3,2:end)') + ...
      bash(7)*func(f_hist(2,1),f_hist(2,2:end)') + bash(8)*func(f_hist(1,1),f_hist(1,2:end)') );

% Calculate the corrector using the Adams-Moulton formula
Y = f_hist(8,2:end)' + h*divisor* ...
    ( moul(1)*func(f_hist(8,1)+h,Y) + moul(2)*func(f_hist(8,1),f_hist(8,2:end)') + ...
      moul(3)*func(f_hist(7,1),f_hist(7,2:end)') + moul(4)*func(f_hist(6,1),f_hist(6,2:end)') + ...
      moul(5)*func(f_hist(5,1),f_hist(5,2:end)') + moul(6)*func(f_hist(4,1),f_hist(4,2:end)') + ...
      moul(7)*func(f_hist(3,1),f_hist(3,2:end)') + moul(8)*func(f_hist(2,1),f_hist(2,2:end)') );

