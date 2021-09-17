classdef clSimulationOutput
   properties
       t double = 0     % seconds from simstart
       Y double = 0     % statevector
       dY double = 0    % differential statevector
       
       grav = struct('spt', [0;0;0],...
           'for', [0;0;0],...
           'tor', [0;0;0])
       
       laser = struct('for', [0;0;0],...
           'tor', [0;0;0],...
           'imp', [0;0;0],...
           'mom', [0;0;0],...
           'projarea', 0)
       
       lt int8 = 0 % 0 no laser engagement, 1 pulsed, 2 cw
       st int8 = 1 % station number
       
       sun = struct('for', [0;0;0],...
           'tor', [0;0;0])
       
       moon = struct('for', [0;0;0],...
           'tor', [0;0;0])
       
       plan = struct('for', [0;0;0],...
           'tor', [0;0;0])
       
       srad = struct('for', [0;0;0],...
           'tor', [0;0;0],...
           'projarea', 0)
       
       albedo = struct('for', [0;0;0],...
           'tor', [0;0;0],...
           'projarea', 0)
       
       infrared = struct('for', [0;0;0],...
           'tor', [0;0;0],...
           'projarea', 0)   
       
       mag = struct('for',[0;0;0])
       
       drag = struct('for', [0;0;0],...
           'tor', [0;0;0],...
           'projarea', 0)
       
       rel = struct('for', [0;0;0],...
           'tor', [0;0;0])
       
       energy = [0;0;0;0]; %work-energy balance REF.: Modeling and Simulation of Long-term
                   %Rotational Dynamics of Large Space Debris Sec. 5.3.2
                   %energy(1): the system/s total mechanical energy
                   %energy(2): the mechanical energy associated with the satellite's orbital motion only
                   %energy(3): the mechanical energy associated with the satellite's rotation
                   %energy(4): the work-energy balance of the satellite's rotational motion
                   %energy(5): the work-energy balance
                   
       incAng = 0; %angle of laser-target vector and the velocity of the target
       
       aux1 = 0
       steptime = 0
       simtime = -1       
   end
end