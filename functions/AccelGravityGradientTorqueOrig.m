function [a_gravity_inertial, g_gravity_body] = AccelGravityGradientTorque(p_ecef, p, ...
    LLA, C_i2b, Inertia, C, S, l_max_a, l_max_g, in_grav_a, in_grav_g)

    p = clPropagator.instance();

    % Use highest maximum degree and order considered for acceleration anf torque
    l_max = l_max_a;
    if (~in_grav_a)
        l_max = l_max_g;
    elseif (in_grav_g && (l_max_g > l_max_a))
        l_max = l_max_g;
    end
    
    % Constants
    mu = 3986004.418E8;
    a = 6378136.3;
    
    % Initialize
    a_gravity_inertial(1:3) = 0;
    g_gravity_body(1:3) = 0;
    
    % Spherical geocentric distance, longitude and latitude (declination) (m and rad)
    r = sqrt(p_ecef(1)*p_ecef(1) + p_ecef(2)*p_ecef(2) + p_ecef(3)*p_ecef(3));
    lon = LLA(3);
    lat = asin(p_ecef(3)/sqrt(p_ecef(1)*p_ecef(1)+p_ecef(2)*p_ecef(2)+p_ecef(3)*p_ecef(3)));
    
    % Legendre Polynomials
%     P(1,1) = 1;
%     P(2,1) = sin(lat);
%     P(2,2) = cos(lat);
%     for l = 0:l_max
%         P(l+1,l+2) = 0;
%         P(l+1,l+3) = 0;
%     end
%     for l = 2:l_max+2
%         P(l+1,1) = ((2*l-1)*sin(lat)*P(l,1)-(l-1)*P(l-1,1))/(l);
%         for m = 1:l
%             P(l+1,m+1) = P(l-1,m+1)+(2*l)*cos(lat)*P(l,m);
%         end
%         P(l+1,l+1) = (2*l-1)*cos(lat)*P(l,l);
%     end
%     
%     %%/ ASPHERICAL ACCELERATION CALCULATION
%     if (in_grav_a)
%     
%         dUdr = 0; % Doesn't take into account spherical term. Counted later in Propagation Function
%         dUd0 = 0;
%         dUdl = 0;
%     
%         % Potential partial derivaties
%         for l = 2:l_max_a
%             for m = 0:l
%                 dUdr = dUdr + (a/r)^l*(l+2)*P(l+1,m+1)*(C(l+1,m+1)*cos(m*lon)+S(l+1,m+1)*sin(m*lon));
%                 dUd0 = dUd0 + (a/r)^l*(P(l+1,m+2)-m*tan(lat)*P(l+1,m+1))*(C(l+1,m+1)*cos(m*lon)+S(l+1,m+1)*sin(m*lon));
%                 dUdl = dUdl + (a/r)^l*m*P(l+1,m+1)*(S(l+1,m+1)*cos(m*lon)-C(l+1,m+1)*sin(m*lon));
%             end
%         end
%         dUdr = -dUdr*mu/(r*r);
%         dUd0 = dUd0*mu/r;
%         dUdl = dUdl*mu/r;
%     
%         %%% Acceleration in TEME Frame
%         a_gravity_inertial(1) = (dUdr/r-p(3)/(r*r*sqrt(p(1)*p(1)+p(2)*p(2)))*dUd0)*p(1) - (dUdl/(p(1)*p(1)+p(2)*p(2)))*p(2);
%         a_gravity_inertial(2) = (dUdr/r-p(3)/(r*r*sqrt(p(1)*p(1)+p(2)*p(2)))*dUd0)*p(2) + (dUdl/(p(1)*p(1)+p(2)*p(2)))*p(1);
%         a_gravity_inertial(3) = dUdr/r*p(3) + sqrt(p(1)*p(1)+p(2)*p(2))/(r*r)*dUd0;
        
        d = norm(r_bf);                     % distance
        latgc = asin(r_bf(3)/d);
        lon = atan2(r_bf(2),r_bf(1));

        [pnm, dpnm] = Legendre(p.AuxParam.n,p.AuxParam.m,latgc);

        dUdr = 0;
        dUdlatgc = 0;
        dUdlon = 0;
        q3 = 0; q2 = q3; q1 = q2;
        for n=0:p.AuxParam.n
            b1 = (-gm/d^2)*(r_ref/d)^n*(n+1);
            b2 =  (gm/d)*(r_ref/d)^n;
            b3 =  (gm/d)*(r_ref/d)^n;
            for m=0:p.AuxParam.m
                q1 = q1 + pnm(n+1,m+1)*(C(n+1,m+1)*cos(m*lon)+S(n+1,m+1)*sin(m*lon));
                q2 = q2 + dpnm(n+1,m+1)*(C(n+1,m+1)*cos(m*lon)+S(n+1,m+1)*sin(m*lon));
                q3 = q3 + m*pnm(n+1,m+1)*(S(n+1,m+1)*cos(m*lon)-C(n+1,m+1)*sin(m*lon));
            end
            dUdr     = dUdr     + q1*b1;
            dUdlatgc = dUdlatgc + q2*b2;
            dUdlon   = dUdlon   + q3*b3;
            q3 = 0; q2 = q3; q1 = q2;
        end

        % Body-fixed acceleration
        r2xy = r_bf(1)^2+r_bf(2)^2;

        ax = (1/d*dUdr-r_bf(3)/(d^2*sqrt(r2xy))*dUdlatgc)*r_bf(1)-(1/r2xy*dUdlon)*r_bf(2);
        ay = (1/d*dUdr-r_bf(3)/(d^2*sqrt(r2xy))*dUdlatgc)*r_bf(2)+(1/r2xy*dUdlon)*r_bf(1);
        az =  1/d*dUdr*r_bf(3)+sqrt(r2xy)/d^2*dUdlatgc;

        a_bf = [ax; ay; az];
    end
    
    
    %%% TORQUE CALCULATION
    if (in_grav_g)

        % Matrices
    
        i2j2 = p(1)*p(1)+p(2)*p(2);
    
        d2rdr2(1:3,1:3) = [...
            p(1)*p(1)-r*r, p(1)*p(2),      p(1)*p(3);...
            p(1)*p(2),     p(2)*p(2)-r*r,  p(2)*p(3);...
            p(1)*p(3),     p(2)*p(3),      p(3)*p(3)-r*r];
        for i=1:3
            for j=1:3
                d2rdr2(i,j) = -d2rdr2(i,j)/(r*r*r);
            end
        end
    
        d20dr2(1:3,1:3) = [...
            p(3)*(2*p(1)^4+p(1)*p(1)*p(2)*p(2)-p(2)^4-p(2)*p(2)*p(3)*p(3)),...
                p(1)*p(2)*p(3)*(3*i2j2+p(3)*p(3)),...
                -p(1)*i2j2*(i2j2-p(3)*p(3));...
            p(1)*p(2)*p(3)*(3*i2j2+p(3)*p(3)),...
                p(3)*(-p(1)^4+p(1)*p(1)*p(2)*p(2)-p(1)*p(1)*p(3)*p(3)+2*p(2)^4),...
                -p(2)*i2j2*(i2j2-p(3)*p(3));...
            -p(1)*i2j2*(i2j2-p(3)*p(3)),...
                -p(2)*i2j2*(i2j2-p(3)*p(3)),...
                -2*p(3)*i2j2*i2j2];
            
        for i=1:3
            for j=1:3
                d20dr2(i,j) = d20dr2(i,j)/(r^4*i2j2^(3/2.0));
            end
        end
    
        d2ldr2(1:3,1:3) = [...
            2*p(1)*p(2),           p(2)*p(2)-p(1)*p(1),    0;...
            p(2)*p(2)-p(1)*p(1),   -2*p(1)*p(2),           0;...
            0,                     0,                      0;];
        
        for i=1:3
            for j=1:3
                d2ldr2(i,j) = d2ldr2(i,j)/(i2j2*i2j2);
            end
        end
    
        drdrT(1:3,1:3) = [ ...
            p(1)*p(1), p(1)*p(2),  p(1)*p(3); ...
            p(1)*p(2), p(2)*p(2),  p(2)*p(3); ...
            p(1)*p(3), p(3)*p(2),  p(3)*p(3)];
        
        for i=1:3
            for j=1:3
                drdrT(i,j) = drdrT(i,j)/(r*r);
            end
        end
    
        d0d0T(1:3,1:3) = [ ...
            p(1)*p(1)*p(3)*p(3),   p(1)*p(2)*p(3)*p(3),    -p(1)*p(3)*i2j2;
            p(1)*p(2)*p(3)*p(3),   p(2)*p(2)*p(3)*p(3),    -p(2)*p(3)*i2j2;
            -p(1)*p(3)*i2j2,       -p(2)*p(3)*i2j2,        i2j2*i2j2];
        
        for i=1:3
            for j=1:3
                d0d0T(i,j) = d0d0T(i,j)/(r*r*r*r*i2j2*i2j2);
            end
        end
    
        dldlT(1:3,1:3) = [ ... 
            p(2)*p(2),     -p(1)*p(2), 0; ...
            -p(1)*p(2),    p(1)*p(1),  0; ...
            0,             0,          0];
        
        for i=1:3
            for j=1:3
                dldlT(i,j) = dldlT(i,j)/(i2j2*i2j2);
            end
        end
    
        drd0T2(1:3,1:3) = [ ...
            -2*p(1)*p(1)*p(3),         -2*p(1)*p(2)*p(3),          p(1)*i2j2-p(1)*p(3)*p(3); ...
            -2*p(1)*p(2)*p(3),         -2*p(2)*p(2)*p(3),          p(2)*i2j2-p(2)*p(3)*p(3); ...
            p(1)*i2j2-p(1)*p(3)*p(3),  p(2)*i2j2-p(2)*p(3)*p(3),   2*p(3)*i2j2 ];
        
        for i=1:3
            for j=1:3
                drd0T2(i,j) = drd0T2(i,j)/(r*r*r*sqrt(i2j2));
            end
        end
    
        drdlT2(1:3,1:3) = [ ...
            -2*p(1)*p(2),          p(1)*p(1)-p(2)*p(2),    -p(2)*p(3); ...
            p(1)*p(1)-p(2)*p(2),   2*p(1)*p(2),            p(1)*p(3); ...
            -p(2)*p(3),            p(1)*p(3),              0];
        
        for i=1:3
            for j=1:3
                drdlT2(i,j) = drdlT2(i,j)/(r*i2j2);
            end
        end
    
        dld0T2(1:3,1:3) = [ ...
            2*p(1)*p(2)*p(3),              p(3)*(p(2)*p(2)-p(1)*p(1)), -p(2)*i2j2; ...
            p(3)*(p(2)*p(2)-p(1)*p(1)),    -2*p(1)*p(2)*p(3),          p(1)*i2j2; ...
            -p(2)*i2j2,                    p(1)*i2j2,                  0];
        
        for i=1:3
            for j=1:3
                dld0T2(i,j) = dld0T2(i,j)/(r*r*i2j2^(3/2.0));
            end
        end
    
        % First and Second partial derivatives
        dUdr = 1;   % Takes into account spherical term
        dUd0 = 0;
        dUdl = 0;
        d2Udr2 = 2; % Takes into account spherical term
        d2Ud02 = 0;
        d2Udl2 = 0;
        d2Ud0dr = 0;
        d2Udldr = 0;
        d2Udld0 = 0;
    
        for l = 2:l_max_g
            for m = 0:l
                dUdr = dUdr + (a/r)^l*(l+2)*P(l+1,m+1)*(C(l+1,m+1)*cos(m*lon)+S(l+1,m+1)*sin(m*lon));
                dUd0 = dUd0 + (a/r)^l*(P(l+1,m+2)-m*tan(lat)*P(l+1,m+1))*(C(l+1,m+1)*cos(m*lon)+S(l+1,m+1)*sin(m*lon));
                dUdl = dUdl + (a/r)^l*m*P(l+1,m+1)*(S(l+1,m+1)*cos(m*lon)-C(l+1,m+1)*sin(m*lon));
            
                d2Udr2 = d2Udr2 + (a/r)^l*(l+2)*(l+1+2)*P(l+1,m+1)*(C(l+1,m+1)*cos(m*lon)+S(l+1,m+1)*sin(m*lon));
                d2Ud02 = d2Ud02 + (a/r)^l*(P(l+1,m+1+2)-(2*m+1)*tan(lat)*P(l+1,m+2)+m*(m*tan(lat)-(1/(cos(lat)*cos(lat))))*P(l+1,m+1))*(C(l+1,m+1)*cos(m*lon)+S(l+1,m+1)*sin(m*lon));
                d2Udl2 = d2Udl2 + (a/r)^l*m*m*P(l+1,m+1)*(C(l+1,m+1)*cos(m*lon)+S(l+1,m+1)*sin(m*lon));
            
                d2Ud0dr = d2Ud0dr + (a/r)^l*(l+2)*(P(l+1,m+2)-m*tan(lat)*P(l+1,m+1))*(C(l+1,m+1)*cos(m*lon)+S(l+1,m+1)*sin(m*lon));
                d2Udldr = d2Udldr + (a/r)^l*m*(l+2)*P(l+1,m+1)*(S(l+1,m+1)*cos(m*lon)-C(l+1,m+1)*sin(m*lon));
                d2Udld0 = d2Udld0 + (a/r)^l*m*(P(l+1,m+2)-m*tan(lat)*P(l+1,m+1))*(S(l+1,m+1)*cos(m*lon)-C(l+1,m+1)*sin(m*lon));
            end
        end
    
        dUdr = -dUdr*mu/(r*r);
        dUd0 = dUd0*mu/r;
        dUdl = dUdl*mu/r;
    
        d2Udr2 = d2Udr2*mu/(r*r*r);
        d2Ud02 = d2Ud02*mu/r;
        d2Udl2 = -d2Udl2*mu/r;
    
        d2Ud0dr = -d2Ud0dr*mu/(r*r);
        d2Udldr = -d2Udldr*mu/(r*r);
        d2Udld0 = d2Udld0*mu/r;
    
        % Acceleration Derivative in TEME Frame
    
        dadr = zeros(3);
        for i=1:3
            for j=1:3
                dadr(i,j) = d2Udr2*drdrT(i,j) + dUdr*d2rdr2(i,j) + d2Ud02*d0d0T(i,j) + d2Udl2*dldlT(i,j) + d2Ud0dr*drd0T2(i,j) + d2Udldr*drdlT2(i,j) + d2Udld0*dld0T2(i,j) + dUd0*d20dr2(i,j) + dUdl*d2ldr2(i,j);
            end
        end
    
        % Rotate into body-fixed frame
%         matrixmult(C_i2b,dadr,dadr_b);
        dadr_b = C_i2b*dadr;
        
        C_b2i = C_i2b';
        
%         matrixmult(dadr_b,C_b2i,G);
        G = dadr_b*C_b2i;
    
        % Gravity-Gradient Torque in body-fixed frame
        g_gravity_body(1) = G(2,3)*(Inertia(3,3)-Inertia(2,2)) - G(1,3)*Inertia(1,2) + G(1,2)*Inertia(1,3) + Inertia(2,3)*(G(2,2)-G(3,3));
        g_gravity_body(2) = G(1,3)*(Inertia(1,1)-Inertia(3,3)) + G(2,3)*Inertia(1,2) - G(1,2)*Inertia(2,3) + Inertia(1,3)*(G(3,3)-G(1,1));
        g_gravity_body(3) = G(1,2)*(Inertia(2,2)-Inertia(1,1)) - G(2,3)*Inertia(1,3) + G(1,3)*Inertia(2,3) + Inertia(1,2)*(G(1,1)-G(2,2));
    end

end