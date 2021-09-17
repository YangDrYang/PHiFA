function t = interpTransmittance(h, za, wl)
    persistent data_t data_h data_wl data_zang F
    
    if isempty(data_t)
        load('lowtran_transmittance.mat', 'data', 'data_heights', 'data_wl_nm', 'data_za');
        data_t = data;
        data_h = data_heights;
        data_wl = data_wl_nm;
        data_zang = data_za;
        F = scatteredInterpolant(data_h(:), data_zang(:), data_wl(:), data_t(:));
        F.Method = 'nearest';
    end
    
    t = F(h/1000, za*180/pi, wl*10^9);
end