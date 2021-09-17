%--------------------------------------------------------------------------
%
% 
%    Adams-Bashforth-Moulton 8th order Integration
%
% Last modified:   2018/01/27   M. Mahooti
%--------------------------------------------------------------------------
clc
clear
format long g

% constants
GM    = 1;                      % gravitational coefficient
e     = 0.1;                    % eccentricity
Kep   = [1, e ,0 ,0 ,0 ,0]';    % (a,e,i,Omega,omega,M)

% header
fprintf( '\nAdams-Bashforth-Moulton 8th order integration\n\n' );

% Initial state of satellite (x,y,z,vx,vy,vz)
y_0 = State(GM, Kep, 0);

% step-size of integration
h = 0.03; % [s]

% initial values
t_0 = 0;
t_end = 86400; % end time [s]

Steps = t_end/h;

f_hist = zeros(8,7);
f_hist(1,1) = t_0;
f_hist(1,2:7) = y_0;
for i=1:7
    y_0 = RK4(@deriv,t_0,y_0,h);
    t_0 = t_0+h;
    f_hist(i+1,1) = t_0;
    f_hist(i+1,2:7) = y_0;
end

for j=8:Steps
    y = ABM8(@deriv,h,f_hist);
    t_0 = t_0+h;
    f_hist(1:7,:) = f_hist(2:8,:);
    f_hist(8,1) = t_0;
    f_hist(8,2:7) = y;    
end
y_ref = State(GM, Kep, t_end); % Reference solution

fprintf(' Accuracy    Digits\n');
fprintf('%10.3g',norm(y-y_ref));
fprintf('%8.2f\n',-log10(norm(y-y_ref)));

