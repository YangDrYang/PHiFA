function eph = eci2ecef_ephemerides(eph, mjd0)

eph_ = eph;
for i = 1:length(eph)
    eph_(2:7,i) = ECI2ECEF(mjd0+eph(1,i)/86400, eph(2:7,i)')';
end


end