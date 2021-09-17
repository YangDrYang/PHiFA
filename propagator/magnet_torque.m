%
%  eddy_torque.c
%  D-SPOSE

%  Created by Luc Sagnieres on 2017-09-12.
%  Copyright 2018 Luc Sagnieres. All rights reserved.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FUNCTION NAME:        eddy_torque.c
%
% DESCRIPTION:          This function calculates the torque due to eddy-
%                       currents (See Eq. (2.20) in Sagnieres (2018)
%                       Doctoral Thesis)
%
% AUTHOR:               Luc Sagnieres
% DATE:                 September 12, 2017
% VERSION:              1
%
% INPUT:                double B_field_b[3]: magnetic field vector in body frame
%                       double B_field_dot_b[3]: time derivative of magnetic field
%                         vector in inertial frame as seen from orbiting
%                         spacecraft, in body-frame components
%                       double w[3]: angular velocity vector in body frame
%                       double magT[3][3]: magnetic tensor
%
% OUTPUT:               double g_eddy[3]: eddy_current torque

function g_eddy = magnet_torque(B_field_b, B_field_dot_b, w, M)
    
    wxB = cross(w,B_field_b);
    Bdiff = wxB-B_field_dot_b;
    
    g_eddy = cross(M*Bdiff,B_field_b);
    
end