%test script

clear;
%% initialice shape

inputstr = clCompoundSegmentsInput;
inputstr(1).bSTL = true;
inputstr(1).stlFile = 'stlfiles/3UCubeSat.stl';
% inputstr(1).stlFile = string('stlfiles\femur.stl');
% inputstr(1).stlFile = string('stlfiles\glove.stl');

inputstr(1).surfAtr = clSA_wilken2015();
a = [-0.01314 0.01753 0.01904 449.81 632.06 913.51];
t = [15.04 1.19 2.41 23.38];
inputstr(1).surfAtr.init(a, t);
inputstr(1).surfAtr.name = 'testsurf1';
inputstr(1).offset = [0; 0; 0];
inputstr(1).name = 'test1';
inputstr(1).density = 2500;
inputstr(1).scale = 1/1000; %units in stl. 1/100 would be centimeter and so on
inputstr(1).resolution = 50;
inputstr(1).solid = true;
inputstr(1).twosided = false;
inputstr(1).thickness = 0.001;

%% initialice target
tmpinput = clCompoundTargetInput();
cmptarget = clCompoundTarget();
cmptarget.init(inputstr, tmpinput);
[vol, mas, moi, pointcloud] = cmptarget.segments.initHollowInertia(); % O(n^3)

%%
plotSegmentWithPointCloud(cmptarget,pointcloud);
% plot2tikz('glove_pointcloud', 1, 'ticks=none');
%%
plotTarget(cmptarget);
% plot2tikz('glove_facets', 1, 'ticks=none');

az = 0:10:360;
for i = 1:length(az)
    view(az(i), 30);
%     saveas(gcf, sprintf('figures/glove_rotate/rot-%d.png', i));
    saveas(gcf, sprintf('figures/3ucubesat_rotate/rot-%d.png', i));
end

