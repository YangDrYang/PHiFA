    % Load IGRF-12 Gauss Coefficients
    function [G,H] = loadMagCoef()
    
    mag_coef = zeros(195,27);
    fid = fopen('propagator/igrf12coeffs.txt','r');
    
    for i =1:4
        fgetl(fid);
    end

    for i = 1:195
        tmp = fscanf(fid,'%c %d %d %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f');
        mag_coef(i,:) = tmp(2:28);
    end
    fclose(fid);
    
    [G,H] = gaus_coef(mag_coef);