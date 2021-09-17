[~,~,isloaded] = pyversion;
if ~isloaded
    pyversion('C:\ProgramData\Miniconda3\envs\py36\python.exe');
end
alt_km = 0;
zenithangle = 0;
p.model=5;
p.h1=alt_km;
p.angle=zenithangle;
p.wlshort= 200;
p.wllong=30000;
T = lowtran_transmission(p);
wl_nm = xarray2mat(T{'wavelength_nm'});


h_vec = -1:0.5:10;
za_vec = 0:5:90;
len_h = length(h_vec);
len_wl = length(wl_nm);
len_za = length(za_vec);

data = zeros(len_h, len_wl,len_za);
data_heights = zeros(len_h, len_wl,len_za);
data_wl_nm = zeros(len_h, len_wl,len_za);
data_za = zeros(len_h, len_wl,len_za);

for za = 1:len_za
    zenithangle = za_vec(za);
    fprintf('Zenithangle: %d\nAltitude: ', zenithangle);
    for h = 1:len_h
        alt_km = h_vec(h);
        fprintf(' %2.0d', alt_km);

        p.model=5;
        p.h1=alt_km;
        p.angle=zenithangle;
        p.wlshort= 200;
        p.wllong=30000;

        T = lowtran_transmission(p);

        trans = squeeze(xarray2mat(T{'transmission'}));
        wl_nm = xarray2mat(T{'wavelength_nm'});

        data(h,:,za) = trans;
        data_heights(h,1:end,za) = h_vec(h);
        data_wl_nm(h,:,za) = wl_nm;
        data_za(h,1:end,za) = za_vec(za);
    end
    fprintf('\n');
end

save('lowtran_transmittance.mat', 'data', 'data_heights', 'data_wl_nm', 'data_za');
figure;
surf(data_wl_nm(:,:,1), data_heights(:,:,1), data(:,:,1));
figure;
surf(data_wl_nm(:,:,end), data_heights(:,:,end), data(:,:,end));