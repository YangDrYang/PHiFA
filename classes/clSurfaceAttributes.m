classdef clSurfaceAttributes < handle
    properties
        %ablation
        thresholdFluence double = -1 % J / cm^2
        %photon pressure
        alpha double = 0.55 % albedo
        beta double = 0.8 % fraction of reflected light that goes into the specular component
        %Ref: Solar Sailing Technology, Dynamics and Mission Applications 1999
        epsf double = 0.05 %front emissivity
%         epsb double = 0.55 %back emissivity
        epsb double = 0.05 %back emissivity
%         Bf double = 0.79 %the front surface will be non-Lambertian, for solar panel only
        Bf double = 0.55 %the front surface will be non-Lambertian, for solar panel only
        Bb double = 0.55 %the back surface will be non-Lambertian, for solar panel only
        %drag - both not needed anymore
        sigma_n double = 0.8 % normal momentum exchange coef
        sigma_t double = 0.8 % tangential momentum exchange coef
        %miscellaneous
        name string = 'none'
    end
    
    methods
        function obj = clSurfaceAttributes()
            %clSurfaceAttributes Construct an instance of this class
        end
    end
    
    methods (Abstract)
        init(obj,varargin)
        cm = getCouplingCoeffientAblation(obj, targetFluence)
    end
end

