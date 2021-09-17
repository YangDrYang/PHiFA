%
%  vectors2angle.c
%  D-SPOSE
%
%  Created by Luc Sagnieres on 2016-02-18.
%  Copyright © 2018 Luc Sagnieres. All rights reserved.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FUNCTION NAME:        vectors2angle.c
%
% DESCRIPTION:          This function calculates the minimum angle between
%                       two vectors
%
% AUTHOR:               Luc Sagnieres
% DATE:                 February 18, 2016
% VERSION:              1
% AUTHOR:               Yang Yang
% DATE:                 July 11, 2019
% VERSION:              2 C -> Matlab
%
% INPUT:                double a[3]: 3x1 vector
%                       double b[3]: 3x1 vector
%
% OUTPUT:               double alpha: minimum angle between two vectors (rad)
%
% COUPLING:             None
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function alpha = vectors2angle(a,b)

    alpha = atan2(norm(cross(a,b)), dot(a,b));
    
end