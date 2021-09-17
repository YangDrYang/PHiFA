function h = getTrueHeight(z,za,h0)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here
pt=[sin(za)*z 6378000+h0+cos(za)*z];
h=norm(pt)-6378000;
end

