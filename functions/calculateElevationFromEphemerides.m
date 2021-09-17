function [elevation, distance] = calculateElevationFromEphemerides(laserobj,tareph)

elevation = zeros(2,length(tareph));
elevation(1,:) = tareph(1,:);

distance = zeros(2,length(tareph));
distance(1,:) = tareph(1,:);

for i = 1:length(tareph)
    laspos = laserobj.getPosition(tareph(1,i));
    las2tar = tareph(2:4,i) - laspos;
%     dotp = dot(las2tar, laspos);
%     ang = acos(dotp/norm(las2tar)/norm(laspos));
    ang = atan2(norm(cross(las2tar,laspos)),dot(las2tar,laspos));
    if ang>pi/2
        elevation(2,i) = 0;
    else
        elevation(2,i) = pi/2 - ang;
    end
    distance(2,i) = norm(las2tar);
end

end
