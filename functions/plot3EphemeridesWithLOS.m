function plot3EphemeridesWithLOS(eph_las,eph_tar)

figure;

% E = wgs84Ellipsoid('kilometres');
% [x,y,z] = ellipsoid(0,0,0,E.SemimajorAxis, E.SemimajorAxis, E.SemiminorAxis);
% surf(x,y,z, 'FaceAlpha', 0.2, 'FaceLighting', 'gouraud', 'EdgeColor', 'none');
% grid;
% hold on;
% 
% state = eph_las(2:4,:)./1000;
% plot3(state(1,:),state(2,:),state(3,:));
% state = eph_tar(2:4,:)./1000;
% plot3(state(1,:),state(2,:),state(3,:));

line = zeros(3,2);
for i = 1:length(eph_las)
    line(:,1) = eph_las(2:4,i);
    dif = abs(eph_tar(1,:)-eph_las(1,i));
    [~, ind] = min(dif);
    r_las2obj = eph_tar(2:4,ind) - eph_las(2:4,i);
%     ang = acos(dot(r_las2obj,eph_las(2:4,i))/norm(r_las2obj)/norm(eph_las(2:4,i)))-pi/2;
    ang = dot(r_las2obj,eph_las(2:4,i));

    if  ang > 0
        line(:,2) = eph_tar(2:4,ind);
        plot3(line(1,:),line(2,:),line(3,:),'b');
        hold on;
    end
end

E = wgs84Ellipsoid('kilometres');
[x,y,z] = ellipsoid(0,0,0,E.SemimajorAxis, E.SemimajorAxis, E.SemiminorAxis);
surf(x,y,z, 'FaceAlpha', 0.2, 'FaceLighting', 'gouraud', 'EdgeColor', 'none');
grid;
hold on;

state = eph_las(2:4,:)./1000;
plot3(state(1,:),state(2,:),state(3,:));
state = eph_tar(2:4,:)./1000;
plot3(state(1,:),state(2,:),state(3,:));

hold off
title('Ephemerides');
xlabel('x [km]');
ylabel('y [km]');
zlabel('z [km]');
axis equal;

end