classdef clCompoundTargetInput
    properties
        % true if inertia calculation can be skipped and the provided data should be used
        bUseProvidedParam logical = false
        inertia(3,3) double = zeros(3) % kg m^2
        mass double        % kg
        barycenter(3,1) double % m
        magnetic_tensor(3,3) double %S*m^4
        
        area_solar double = 110.5 % effective solar area
        area_drag double = 62.5 % effective drag area
        Cr double = 1.0 % solar radiation pressure coefficient
        Cd double = 2.0 % drag coefficient
    end
end