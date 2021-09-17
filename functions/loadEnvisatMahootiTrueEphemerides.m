function eph = loadEnvisatMahootiTrueEphemerides()
    run propagator/True_EnvisatStates.m
    [cTE, rTE] = size(True_Eph);
    eph = zeros(cTE,rTE+1);
    eph(:,1) = 0:60:60*(cTE-1);
    eph(:,2:7) = True_Eph;
    eph = eph';
end