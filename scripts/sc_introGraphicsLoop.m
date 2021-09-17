load logfiles/sim0151.mat

t_step = 10*30;
t_interval = 2*t_step;

[pert, prop, vars] = reformatSimulationData(initstate, propagator.output_table);
steph = propagator.station.ephemerides;
steph(2:7,:) = steph(2:7,:)./1000;


t = 0:t_step:49000;

E = wgs84Ellipsoid('kilometres');
[x,y,z] = ellipsoid(0,0,0,E.SemimajorAxis, E.SemimajorAxis, E.SemiminorAxis);

prop(2:7,:) = prop(2:7,:)./1000;

for i = 1:length(t)
    t_start = t(i)-t_interval/2;
    t_end = t(i)+t_interval/2;
    
    if t_start<prop(1,1)
        t_start = prop(1,1);
    end
    
    if t_end>prop(1,end)
        t_end = prop(1,end);
    end
    
    i_dot = find(prop(1,:)==t(i));
    i_start = find(prop(1,:)==t_start);
    i_end = find(prop(1,:)==t_end);
    il_dot = find(steph(1,:)==t(i));
    il_start = find(steph(1,:)==t_start);
    il_end = find(steph(1,:)==t_end);
    
    figure;
    surf(x,y,z, 'FaceAlpha', 0.2, 'FaceLighting', 'gouraud', 'EdgeColor', 'none');
    hold on;
    
    plot3(prop(2,i_start:i_end),prop(3,i_start:i_end),prop(4,i_start:i_end));
    plot3(prop(2,i_dot),prop(3,i_dot),prop(4,i_dot),'k.');
    plot3(steph(2,il_start:il_end),steph(3,il_start:il_end),steph(4,il_start:il_end));
    plot3(steph(2,il_dot),steph(3,il_dot),steph(4,il_dot),'k.');
    
    if pert(2, i_dot)~=0
        line = [prop(2:4,i_dot), steph(2:4,il_dot)];
        plot3(line(1,:),line(2,:),line(3,:));
    end
    
    view([310 45]);
    axis equal;
    xlabel('x [km]');
    ylabel('y [km]');
    zlabel('z [km]');
    saveas(gcf,sprintf('figures/intro/img-%d.png',i));
    close gcf
end