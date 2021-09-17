classdef clSA_const < clSurfaceAttributes
    %clSA_const Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        const_cm
    end
    
    methods
        function obj = clSA_const()
            %clSA_const Constant coupling coefficient
        end
        
        function init(obj,ccm)
            obj.const_cm = ccm;
        end
        
        function cm = getCouplingCoeffientAblation(obj, ~)
            cm = obj.const_cm;
        end
    end
end

