function plot3Ephemerides(eph,station)

figure;

E = wgs84Ellipsoid('kilometres');
[x,y,z] = ellipsoid(0,0,0,E.SemimajorAxis, E.SemimajorAxis, E.SemiminorAxis);
h1 = surf(x,y,z, 'FaceAlpha', 0.2, 'FaceLighting', 'gouraud', 'EdgeColor', 'none');
legend('Earth');
hold on;

% if initstate.sec > 0
%     mjd0 = Mjday(initstate.year,initstate.mon,initstate.day,...
%         initstate.hour,initstate.min,initstate.sec);
% else
%     mjd0 = initstate.mjd;
% end
% 
% eph_ecef = eci2ecef_ephemerides(eph, mjd0);
% state = eph_ecef(2:4,:)./1000;
state = eph(2:4,:)./1000;
% [h1,h2] = EarthOrbitPlot(state);
h2 = plot3(state(1,1:end),state(2,1:end),state(3,1:end));
hold on
plot3(state(1,:),state(2,:),state(3,:));
h3 = plot3(state(1,1),state(2,1),state(3,1),'k*');
plot3(state(1,1000),state(2,1000),state(3,1000),'k>');
h4 = plot3(state(1,end),state(2,end),state(3,end),'ko');
% plot3(state(1,1000),state(2,1000),state(3,1000),'k+');
% plot3(state(1,2000),state(2,2000),state(3,2000),'k+');

if ~station.bGroundBased
    ind = 31:7231;
%     eph_ecef = eci2ecef_ephemerides(station.ephemerides, mjd0);
%     state = eph_ecef(2:4,ind)./1000;
    state = station.ephemerides(2:4,ind)./1000;
    h5 = plot3(state(1,:),state(2,:),state(3,:));
    h6 = plot3(state(1,1),state(2,1),state(3,1),'r*');
    plot3(state(1,1000),state(2,1000),state(3,1000),'r>');
%     plot3(state(1,1000),state(2,1000),state(3,1000),'r+');
%     plot3(state(1,2000),state(2,2000),state(3,2000),'r+');
    h7 = plot3(state(1,end),state(2,end),state(3,end),'ro');
end

hold off;
% title('Ephemerides');
xlabel('x [km]','FontSize',14);
ylabel('y [km]','FontSize',14);
zlabel('z [km]','FontSize',14);
legend([h1,h2,h3,h4,h5,h6,h7],'Earth','Target Object (TO)','TO Staring Point','TO Ending Point','Laser Station (LS)','LS Starting Point',...
    'LS Ending Point','FontSize',14)
% legend([h2,h3,h4,h5,h6,h7],'Laser Station (LS)','LS Starting Point',...
%     'LS Ending Point','Target Object (TO)','TO Staring Point','TO Ending Point','FontSize',14)
axis equal;
grid on;
end