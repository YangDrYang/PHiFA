function evalHitMethods(out_beam, out_net, out_area, res, is)

tolerance = 0.05;
spacing(1,length(res)) = 0;
imp_beam(4,length(out_beam)) = 0;
for_beam(4,length(out_beam)) = 0;
imp_net(4,length(out_beam)) = 0;
for_net(4,length(out_beam)) = 0;
imp_area(4,length(out_beam)) = 0;
for_area(4,length(out_beam)) = 0;
time_taken_beam(1,length(out_beam)) = 0;
time_taken_net(1,length(out_beam)) = 0;
time_taken_area(1,1:length(out_beam)) = out_area.aux1;

for i = 1:length(out_beam)
    spacing(i) = 1/res(i);
    if is.a_lase
        imp_beam(1:3,i) = out_beam(i).laser.imp;
        imp_beam(4,i) = norm(out_beam(i).laser.imp);
        for_beam(1:3,i) = out_beam(i).laser.for;
        for_beam(4,i) = norm(out_beam(i).laser.for);
        
        imp_net(1:3,i) = out_net(i).laser.imp;
        imp_net(4,i) = norm(out_net(i).laser.imp);
        for_net(1:3,i) = out_net(i).laser.for;
        for_net(4,i) = norm(out_net(i).laser.for);
        
        for_area(1:3,i) = out_area.laser.for;
        imp_area(1:3,i) = out_area.laser.imp;        
        for_area(4,i) = norm(out_area.laser.for);
        imp_area(4,i) = norm(out_area.laser.imp);
    end
    
    time_taken_beam(i) = out_beam(i).aux1;
    time_taken_net(i) = out_net(i).aux1;

end

if is.a_lase
    figure;
    semilogx(spacing, imp_beam(4,:));
    hold on;
    semilogx(spacing, imp_net(4,:));
    semilogx(spacing, imp_area(4,:));
    
    yline(imp_beam(4,end)+imp_beam(4,end)*tolerance,'HandleVisibility','off');
    yline(imp_beam(4,end)-imp_beam(4,end)*tolerance,'HandleVisibility','off');
    yline(imp_net(4,end)+imp_net(4,end)*tolerance,'HandleVisibility','off');
    yline(imp_net(4,end)-imp_net(4,end)*tolerance,'HandleVisibility','off');
    yline(imp_area(4,end)+imp_area(4,end)*tolerance,'HandleVisibility','off');
    yline(imp_area(4,end)-imp_area(4,end)*tolerance,'HandleVisibility','off');
    
    xlabel('Resolution [m]');
    ylabel('Laser Impulse [Ns]');
    legend('beam','net','area');
    title('Laser Impulse');
    grid;
    hold off;

    figure;
    semilogx(spacing, for_beam(4,:));
    hold on;
    semilogx(spacing, for_net(4,:));
    semilogx(spacing, for_area(4,:));
    
    yline(for_beam(4,end)+for_beam(4,end)*tolerance,'HandleVisibility','off');
    yline(for_beam(4,end)-for_beam(4,end)*tolerance,'HandleVisibility','off');
    yline(for_net(4,end)+for_net(4,end)*tolerance,'HandleVisibility','off');
    yline(for_net(4,end)-for_net(4,end)*tolerance,'HandleVisibility','off');
    yline(for_area(4,end)+for_area(4,end)*tolerance,'HandleVisibility','off');
    yline(for_area(4,end)-for_area(4,end)*tolerance,'HandleVisibility','off');
    
    xlabel('Resolution [m]');
    ylabel('Laser Force (CW) [N]');
    legend('beam','net','area');
    title('CW Force');
    grid;
    hold off;
end

figure;
semilogx(spacing, time_taken_beam);
hold on;
semilogx(spacing, time_taken_net);
semilogx(spacing, time_taken_area);
xlabel('Resolution [m]');
ylabel('Time taken [s]');
title('Time consumption');
legend('beam','net','area');
grid;
hold off;

end