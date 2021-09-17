classdef clRay < handle
    %clRay Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        origin(3,1) double = [0; 0; 0] % m
        direction(3,1) double = [1; 0; 0] % m
        
        xdir(3,1) double = [0; 1; 0]
        ydir(3,1) double = [0; 0; 1]
    end
end

