% sc_plotTumod.m

figure;
tumod_d = @tumodSLCDay;
x = 0:50:30000;
y_d = tumod_d(x,0,0);
semilogy(x,y_d);
grid;
hold on
tumod_n = @tumodSLCNight;
yn = tumod_n(x,0,0);
semilogy(x,yn);
yp4 = tumod_d(x,pi/4,0);
semilogy(x,yp4);
yp289 = tumod_d(x,pi/2*80/90,0);
semilogy(x,yp289);
xlim([0 3E4]);
ylim([1E-18 1E-13]);
xlabel('Distance from laser station [m]');
ylabel('Turbulence-Parameter C_n^2 [-]');
legend('SLC Day, zenith angle 0°', ...
    'SLC Night, zenith angle 0°', ...
    'SLC Day, zenith angle 45°', ...
    'SLC Day, zenith angle 80°');

plot2tikz('cn2', 0.6);