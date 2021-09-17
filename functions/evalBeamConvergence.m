function evalBeamConvergence(out, res, is, figname)

laser_imp(4,length(out)) = 0;
laser_for(4,length(out)) = 0;
drag(4,length(out)) = 0;
srad(4,length(out)) = 0;
time_taken(1,length(out)) = 0;

for i = 1:length(out)
    if is.a_lase
        laser_imp(1:3,i) = out(i).laser.imp;
        laser_imp(4,i) = norm(out(i).laser.imp);
        laser_for(1:3,i) = out(i).laser.for;
        laser_for(4,i) = norm(out(i).laser.for);
    end
    if is.a_drag
        drag(1:3,i) = out(i).drag.for;
        drag(4,i) = norm(out(i).drag.for);
    end
    if is.a_srad
        srad(1:3,i) = out(i).srad.for;
        srad(4,i) = norm(out(i).srad.for);
    end
    
    time_taken(i) = out(i).aux1;
end

if is.a_lase
    figure('Name', figname);
    plot(res, laser_imp(1,:));
    hold on;
    plot(res, laser_imp(2,:));
    plot(res, laser_imp(3,:));
    plot(res, laser_imp(4,:));
    xlabel('Beam resolution [beams/m]');
    ylabel('Laser Impulse [Ns]');
    legend('x','y','z','norm');
    title('Laser Impulse');
    grid;
    hold off;

    figure('Name', figname);
    plot(res, laser_for(1,:));
    hold on;
    plot(res, laser_for(2,:));
    plot(res, laser_for(3,:));
    plot(res, laser_for(4,:));
    xlabel('Beam resolution [beams/m]');
    ylabel('Laser Force (CW) [N]');
    legend('x','y','z','norm');
    title('CW Force');
    grid;
    hold off;
end

if is.a_drag
    figure('Name', figname);
    plot(res, drag(1,:));
    hold on;
    plot(res, drag(2,:));
    plot(res, drag(3,:));
    plot(res, drag(4,:));
    xlabel('Beam resolution [beams/m]');
    ylabel('Drag [N]');
    legend('x','y','z','norm');
    title('Drag');
    grid;
    hold off;
end

if is.a_srad
    figure('Name', figname);
    plot(res, srad(1,:));
    hold on;
    plot(res, srad(2,:));
    plot(res, srad(3,:));
    plot(res, srad(4,:));
    xlabel('Beam resolution [beams/m]');
    ylabel('SRP [N]');
    legend('x','y','z','norm');
    title('Solar Radiation Pressure');
    grid;
    hold off;
end

figure('Name', figname);
plot(res, time_taken);
xlabel('Beam resolution [beams/m]');
ylabel('Time taken for all [s]');
title('Time consumption');
grid;

end