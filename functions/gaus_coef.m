%
%  gaus_coef.c
%  D-SPOSE
%
%  Created by Luc Sagnieres on 2017-08-22.
%  Copyright © 2018 Luc Sagnieres. All rights reserved.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FUNCTION NAME:        gaus_coef.c
%
% DESCRIPTION:          This function will reformat the IGRF-12 coefficients
%
% AUTHOR:               Luc Sagnieres
% DATE:                 August 22, 2017
% VERSION:              1
% AUTHOR:               Yang Yang
% DATE:                 July 09, 2019
% VERSION:              2 C-> Matlab
%
% INPUT:                double mag_coef[195][27]: IGRF-12 coefficients form input
%
% OUTPUT:               double G[14][14][25]: magnetic potential coefficients
%                       double H[14][14][25]: magnetic potential coefficients
%
% COUPLING:             None
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [G,H] = gaus_coef(mag_coef)
G = zeros(14,14,25);
H = zeros(14,14,25);
k = 1;

for i = 1:13
    for j = 1:25
        G(mag_coef(k,1)+1,mag_coef(k,2)+1,j) = mag_coef(k,j+2);
    end
    k = k+1;
    for l = 1:i
        for j = 1:25
            G(mag_coef(k,1)+1,mag_coef(k,2)+1,j) = mag_coef(k,j+2);
        end
        k = k+1;
        for j = 1:25
            H(mag_coef(k,1)+1,mag_coef(k,2)+1,j) = mag_coef(k,j+2);
        end
        k = k+1;
    end
end
            
