classdef clSA_wilken2015 < clSurfaceAttributes
    %clSA_wilken2015 Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        %ablation
        a(1,6) double = [0 0 0 0 0 0]
        t(1,4) double = [0 0 0 0]
    end
    
    methods
        function obj = clSA_wilken2015()
            %clSA_wilken2015 Construct an instance of this class
            %   Detailed explanation goes here
        end
        
        function init(obj,a,t)
            obj.a = a;
            obj.t = t;
            obj.thresholdFluence = - a(3) * log( 1 + a(1)/a(2) );
        end
        
        function cm = getCouplingCoeffientAblation(obj, targetFluence)
            if targetFluence>obj.thresholdFluence
                massLossPerSurfacearea = obj.a(1) + obj.a(2)*(1 - exp(-targetFluence/obj.t(1))) + obj.a(3)*(1 - exp(-targetFluence/obj.t(1)));
                specImpulse = obj.a(4) + obj.a(5)*(1 - exp(-targetFluence/obj.t(3))) + obj.a(5)*(1 - exp(-targetFluence/obj.t(4)));
                cm = (massLossPerSurfacearea * 9.81 * specImpulse)/targetFluence;
            else
                cm = 0;
            end
        end
    end
end

