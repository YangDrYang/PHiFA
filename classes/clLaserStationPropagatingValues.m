classdef clLaserStationPropagatingValues
    properties
        % true if laserstations pass is to propagate prior to target and
        % not parallel
        pre_propagate = true
        
        area_solar % effective solar area
        area_drag % effective drag area
        mass % mass
        Cr % solar radiation pressure coefficient
        Cd % drag coefficient
    end
end