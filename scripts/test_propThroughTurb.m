% for i=12:18
load test_propThroughTurb.mat
%% init
last = clLaserStation();

last.referencePower = 10e3;
last.Aperture = 1;
last.FWHM = 0.5;
last.beamType = eLaserBeam.Gaussian;
last.tumod = @tumodSLCDay;
last.lla(3)=0;

last.xv = [0 6371*sin(pi/4) 6371*cos(pi/4) 0 0 0];

last.init();

%% target

if ~exist('cmptarget', 'var')
    load('var_cmptarget.mat');
end

cmptarget.xv = [0 0 400000 7000 0 0];

%% calc distrib
tic
last.setNextPulse(cmptarget);
%% plotting distrib
% close all;
figure;
subplot(1,3,1);
magnpt = abs(last.Uturb);
contour(last.xn,last.yn,magnpt);
shading interp;
h = colorbar;
xlabel('x [m]');
ylabel('y [m]');
ylabel(h, 'Intensity [W/m^2]');
title('Propagation through Turbulence');
subplot(1,3,2);
magnpt = abs(last.Uvac);
contour(last.xn,last.yn,magnpt);
shading interp;
h = colorbar;
xlabel('x [m]');
ylabel('y [m]');
ylabel(h, 'Intensity [W/m^2]');
title('Propagation through Vacuum');
subplot(1,3,3);
magnpt = abs(last.Uin);
contour(last.x1,last.y1,magnpt);
shading interp;
h = colorbar;
xlabel('x [m]');
ylabel('y [m]');
ylabel(h, 'Intensity [W/m^2]');
title('Laserbeam leaving Telescope');

savefig('figures\test_propThroughTurb.fig');