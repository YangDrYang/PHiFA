function distance = point_to_plane(planeBase,planeNormal,pt)
%point_to_plane Summary of this function goes here
%   Detailed explanation goes here
distance = abs(dot(planeBase-pt,planeNormal))/norm(planeNormal);
end

