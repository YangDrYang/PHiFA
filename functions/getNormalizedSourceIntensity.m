function n_int = getNormalizedSourceIntensity(shaMod, rSC, rSource)

if shaMod == 1
    if dot(rSC, rSource)<0
        p = clPropagator.instance();

        s_rel = rSource - rSC;
        dist = norm(s_rel);
        %check for shadow
        %assuming perfectly cylindrical shadow cone \citep{frueh2013}
        refdist = sin(acos(dot(s_rel, rSC)/...
            (dist*norm(rSC))))*norm(rSC); %derived from equation 8
        if refdist<p.const.R_Earth  %shadow
            n_int = 0;
        else                    %light
            n_int = 1;
        end
    else
        n_int = 1;
    end
    
%     n_int = Cylindrical(rSC, rSource); % uses Mahooti Function
elseif shaMod == 2
    [n_int, ~] = Eclipse( rSC./1000, rSource./1000 );
elseif shaMod == 3
    n_int = shadow(rSource,rSC);
end

end
