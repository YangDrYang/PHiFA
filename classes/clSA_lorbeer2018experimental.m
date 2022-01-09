classdef clSA_lorbeer2018experimental < clSurfaceAttributes
    %clSA_lorbeer2018experimental Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        % empirical fit
        phi0 double
        b double
        c double
        deltaphi double
        
        % setup
        wavelength double
        pulselength double
        A double
        
        % deviation
        bUseDeviation logical = false
    end
    
    methods
        function obj = clSA_lorbeer2018experimental()
            %clSA_lorbeer2018experimental Construct an instance of this class
            %   Detailed explanation goes here
        end
        
        function init(obj, tau, wl, A)
            obj.wavelength = wl;
            obj.pulselength = tau;
            obj.A = A;
            tau = tau*10^9;
            [obj.phi0,obj.deltaphi,obj.b,obj.c] = ...
                loadSAParameters_lorbeer2018experimental('./inputfiles/lorbeer2018experimental_data.txt', tau);
            obj.thresholdFluence = obj.phi0(1)*10000;
        end
        
        function cm = getCouplingCoeffientAblation(obj, targetFluence)
            if targetFluence>obj.thresholdFluence
                if obj.bUseDeviation
                    phi0_rand = randn*obj.phi0(2);
                    deltaphi_rand = randn*obj.deltaphi(2);
                    b_rand = randn*obj.b(2);
                    c_rand = randn*obj.c(2);
                else
                    phi0_rand = 0;
                    deltaphi_rand = 0;
                    b_rand = 0;
                    c_rand = 0;
                end
                targetFluence = targetFluence/10000; %conversion to J/cm^2
                cm = ((targetFluence - (obj.phi0(1)+phi0_rand))/((obj.deltaphi(1)+deltaphi_rand) + (targetFluence - (obj.phi0(1)+phi0_rand))) ...
                    * (obj.b(1)+b_rand) * 12.46 * obj.A^(7/16) * ...
                    ( sqrt(obj.pulselength)/(obj.wavelength * targetFluence ) )^(obj.c(1)+c_rand) );
                if cm<0
                    cm=0;
                end
            else
                cm = 0;
            end
        end
    end
end

