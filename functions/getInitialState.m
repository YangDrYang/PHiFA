function Y0 = getInitialState(raan, inclination, height)
    raan = raan*pi/180;
    inclination = inclination*pi/180;
    p = clPropagator.instance();
    r = p.const.R_Earth + height*1000;
    v = sqrt(p.const.GM_Earth/r);
    
    Y0 = [  r*cos(raan);r*sin(raan);0; ...
            -v*cos(inclination)*sin(raan);v*cos(raan)*cos(inclination);v*sin(inclination) ];
end