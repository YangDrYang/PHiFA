function plot_transmission(T)
trans = squeeze(xarray2mat(T{'transmission'})).*100;
wl_nm = xarray2mat(T{'wavelength_nm'});

figure
plot(wl_nm, trans(:,1))
hold on
plot(wl_nm, trans(:,2))
plot(wl_nm, trans(:,3))
% ylim([1e-4,1])
xlim([200 3000]);
xlabel('Wavelength [nm]')
ylabel('Transmittance [%]')
hold off
grid;
pbaspect([6 2 2]);
end