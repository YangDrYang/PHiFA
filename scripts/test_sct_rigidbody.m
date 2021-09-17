inertia = [2000 300 -200;...
300 4000 -100;...
-200 -100 1000];
invInertia = inv(inertia);
tEnd = 1000;
dT = 1.0;
nSim = tEnd/dT;
torque = [0;0;0.1];
hPlot = zeros(3,nSim);
x = [1;0;0;0;0;0;0];

for k = 1:nSim
hPlot(:,k) = QTForm(x(1:4),inertia*x(5:7));
x = RK4( 'FRB', x, dT, 0, inertia,...
invInertia, QForm( x(1:4), torque ) );
end
Plot2D( (0:(nSim-1))*dT, hPlot, 'Time (sec)',...
['Hx';'Hy';'Hz'],sprintf('Rigid Body Dynamics for dT = %1.1f s',dT));
PrintFig(1,1,1,'RBSim_1p0');

% dT = 0.2;
% nSim = tEnd/dT;
% hPlot = zeros(3,nSim);
% x = [1;0;0;0;0;0;0];
% 
% for k = 1:nSim
% hPlot(:,k) = QTForm(x(1:4),inertia*x(5:7));
% x = RK4( ’FRB’, x, dT, 0, inertia,
% invInertia, QForm( x(1:4), torque ) );
% end
% Plot2D( (0:(nSim-1))*dT, hPlot, ’Time (sec)’,
% [’Hx’;’Hy’;’Hz’],sprintf(’Rigid Body
% Dynamics for dT = %1.1f s’,dT))
% PrintFig(1,1,2,’RBSim_0p2’)