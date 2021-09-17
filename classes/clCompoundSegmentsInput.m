classdef clCompoundSegmentsInput
    %clCompoundSegmentsInput Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        bSTL logical       % true if input facets come from stl file
        stlFile char       % stl file name
        segobj clCompoundSegment   % segment holding facets
        
        surfAtr clSurfaceAttributes = clSA_wilken2015()     % surface Attributes
        offset(3,1) double % offset to object from initial zero
        name char          % arbituary name
        density double     % density in kg per m^3
        scale double       % scale to bring length unit to a meter. depends on stl file
        resolution double  % resolution in cells per meter
        solid logical      % true of solid else false
        twosided logical   % if elements can be considered twosided true, else normal vector in stl file defines outward - doenst do anything right now
        thickness double   % thickness in m
        
        % true if inertia calculation can be skipped and the provided data should be used
        bUseProvidedParam logical = false
        inertia(3,3) double = zeros(3) % kg m^2
        mass double        % kg
    end
end

