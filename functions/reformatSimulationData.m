function [pert, prop, vars] = reformatSimulationData(is, sim_out)
% function reformatSimulationData Takes in output_table of propagator
% object and creates easy to handle matrices
%
% pert - perturbations
% pert (1,:) seconds after simulation start
% pert(2:4, :) = laser impulse
% pert(5:7, :) = laser angular momentum
% pert(8:10, :) = laser cw force
% pert(11:13, :) = laser cw torque
% pert(14:16, :) = gravitation
% pert(17:19, :) = gravity torque
% pert(20:22, :) = drag
% pert(23:25, :) = drag torque
% pert(26:28, :) = SRP
% pert(29:31, :) = SRP torque
%
% prop - propagation -> ephemerides
% prop(2:4,:) = position
% prop(5:7,:) = velocity
% prop(8:11,:) = quaternion
% prop(12:14,:) = angular velocity
%
% vars - other variables
% vars(2,:) - kind of laser engagement
% 0 - no laser engagement
% 1 - cw
% 2 - pulsed


time(1,length(sim_out)-1) = 0;
laser_imp(3,length(sim_out)-1) = 0;
laser_for(3,length(sim_out)-1) = 0;
grav(3,length(sim_out)-1) = 0;
srad(3,length(sim_out)-1) = 0;
drag(3,length(sim_out)-1) = 0;
laser_mom(3,length(sim_out)-1) = 0;
laser_tor(3,length(sim_out)-1) = 0;
gravt(3,length(sim_out)-1) = 0;
sradt(3,length(sim_out)-1) = 0;
dragt(3,length(sim_out)-1) = 0;
eph(14,length(sim_out)-1) = 0;

vars = zeros(2, length(sim_out)-1);
pert = zeros(31, length(sim_out)-1);

for i = 2:length(sim_out)
    if i~=1 && sim_out(i).t==0
        time = time(:,1:i-2);
        laser_imp = laser_imp(:,1:i-2);
        laser_for = laser_for(:,1:i-2);
        grav = grav(:,1:i-2);
        srad = srad(:,1:i-2);
        drag = drag(:,1:i-2);
        laser_mom = laser_mom(:,1:i-2);
        laser_tor = laser_tor(:,1:i-2);
        gravt = gravt(:,1:i-2);
        sradt = sradt(:,1:i-2);
        dragt = dragt(:,1:i-2);
        eph = eph(:,1:i-2);
        pert = pert(:,1:i-2);
        vars = vars(:,1:i-2);
        break;
    end
    time(i-1) = sim_out(i).t(end);
    eph(1,i-1) = time(i-1);
    
    vars(1,i-1) = time(i-1);
    vars(2,i-1) = sim_out(i).lt;
    
    eph(2:14,i-1) = sim_out(i).Y;
    if is.a_lase
        laser_imp(1:3,i-1) = sim_out(i).laser.imp;
        laser_imp(4,i-1) = norm(sim_out(i).laser.imp);
        laser_for(1:3,i-1) = sim_out(i).laser.for;
        laser_for(4,i-1) = norm(sim_out(i).laser.for);
    else
        laser_imp(1:4,i-1) = [0;0;0;0];
        laser_for(1:4,i-1) = [0;0;0;0];
    end
    if is.a_drag
        drag(1:3,i-1) = sim_out(i).drag.for;
        drag(4,i-1) = norm(sim_out(i).drag.for);
    else
        drag(1:4,i-1) = [0;0;0;0];
    end
    if is.a_srad
        srad(1:3,i-1) = sim_out(i).srad.for;
        srad(4,i-1) = norm(sim_out(i).srad.for);
    else
        srad(1:4,i-1) = [0;0;0;0];
    end
    if is.a_grav
        grav(1:3,i-1) = sim_out(i).grav.for;
        grav(4,i-1) = norm(sim_out(i).grav.for);
    else
        grav(1:4,i-1) = [0;0;0;0];
    end
    
    if is.g_lase
        laser_mom(1:3,i-1) = sim_out(i).laser.mom;
        laser_mom(4,i-1) = norm(sim_out(i).laser.mom);
        laser_tor(1:3,i-1) = sim_out(i).laser.tor;
        laser_tor(4,i-1) = norm(sim_out(i).laser.tor);
    else
        laser_mom(1:4,i-1) = [0;0;0;0];
        laser_tor(1:4,i-1) = [0;0;0;0];
    end
    if is.g_drag
        dragt(1:3,i-1) = sim_out(i).drag.tor;
        dragt(4,i-1) = norm(sim_out(i).drag.tor);
    else
        dragt(1:4,i-1) = [0;0;0;0];
    end
    if is.g_srad
        sradt(1:3,i-1) = sim_out(i).srad.tor;
        sradt(4,i-1) = norm(sim_out(i).srad.tor);
    else
        sradt(1:4,i-1) = [0;0;0;0];
    end
    if is.g_grav
        gravt(1:3,i-1) = sim_out(i).grav.tor;
        gravt(4,i-1) = norm(sim_out(i).grav.tor);
    else
        gravt(1:4,i-1) = [0;0;0;0];
    end
end

pert(1, :) = time;

nrows = 2;
pert(nrows:nrows+2, :) = laser_imp(1:3,:);
pert(nrows+3:nrows+5, :) = laser_mom(1:3,:);
nrows = nrows+6;
pert(nrows:nrows+2, :) = laser_for(1:3,:);
pert(nrows+3:nrows+5, :) = laser_tor(1:3,:);
nrows = nrows+6;
pert(nrows:nrows+2, :) = grav(1:3,:);
nrows = nrows+3;
pert(nrows:nrows+2, :) = gravt(1:3,:);
nrows = nrows+3;
pert(nrows:nrows+2, :) = drag(1:3,:);
nrows = nrows+3;
pert(nrows:nrows+2, :) = dragt(1:3,:);
nrows = nrows+3;
pert(nrows:nrows+2, :) = srad(1:3,:);
nrows = nrows+3;
pert(nrows:nrows+2, :) = sradt(1:3,:);

% prop = zeros(size(eph,1)+1,length(eph));
prop = eph;

end