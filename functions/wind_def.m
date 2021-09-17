function v_rel_i = wind_def(mjd, xv_ecef)
%% function calulates relative wind vector in inertial frame given position vector of 
% object in inertial frame and MJD

    %% insert WIND MODEL here
    wind_hor_ecef = [0;0;0];
    
    %% Co-rotating Winds -> transforming to inertial
    tmp = [xv_ecef(1:3); wind_hor_ecef];
    winds = ECEF2ECI(mjd, tmp');
    v_rel_i = winds(4:6);

end