function eph = ecef2eci_ephemerides(eph, mjd0)

for i = 1:length(eph)
    eph(2:7,i) = ECEF2ECI(mjd0+eph(1,i)/86400, eph(2:7,i)')';
end


end