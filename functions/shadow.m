function n_int = shadow(r_sun, p)

n_int = 1;

alpha_umb = 0.264121687*pi/180.0;
alpha_pen = 0.269007205*pi/180.0;

rrd = dotprod(p,r_sun);
rn = norm(p);

if rrd < 0
    r_sun_neg = -r_sun;

    sigma = atan2d(norm(cross(r_sun_neg,p)),dot(r_sun_neg,p));
    sat_hori = rn*cos(sigma);
    sat_vert = rn*sin(sigma);
    R_earth = 6378.137e3;

    xx = R_earth/sin(alpha_pen);

    pen_vert = tan(alpha_pen)*(xx+sat_hori);

    if sat_vert<=pen_vert

        yy = R_earth/sin(alpha_umb);
        umb_vert = tan(alpha_umb)*(yy-sat_hori);

        if sat_vert<=umb_vert
            n_int = 0;
        else

            R_sun = 6.957e8;
            rsat = norm(p);

            rsr = r_sun - p;
            nrsr = norm(rsr);

            a = asin(R_sun/nrsr);

            b = asin(R_earth/rsat);

            sdotr = dotprod(p,rsr);

            c = acos(-sdotr/(rsat*nrsr));

            x = (c*c + a*a - b*b)/(2*c);
            y = sqrt(a*a-x*x);

            A = a*a*acos(x/a)+b*b*acos((c-x)/b)-c*y;

            if ((a+b) <= c)
                n_int = 1;
            elseif (c<(b-a))||(c<(a-b))
                n_int = 0;
            else
                n_int = 1-A/(pi*a*a);
            end
            
        end
        
    end
end
    
end