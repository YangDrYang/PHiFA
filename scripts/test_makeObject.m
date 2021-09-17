clear;
%% initialice shape

inputstr = clCompoundSegmentsInput;

inputstr(1).bSTL = false;
inputstr(1).stlFile = "";

inputstr(1).surfAtr = clSA_wilken2015();
a = [-0.01314 0.01753 0.01904 449.81 632.06 913.51];
t = [15.04 1.19 2.41 23.38];
inputstr(1).surfAtr.init(a, t);
inputstr(1).surfAtr.name = 'testsurf1';
inputstr(1).offset = [0 0 0];
inputstr(1).name = 'test1';
inputstr(1).density = 2500;
inputstr(1).scale = 1/1000; %units in stl. 1/100 would be centimeter and so on
inputstr(1).resolution = 200;
inputstr(1).solid = true;
inputstr(1).twosided = false;
inputstr(1).thickness = 0.001;

%% initialice target
cmptarget = clCompoundTarget();

%% object
% testseg = makeCube(0.1);
% testseg = makePlate(0.25, 0.1);
% testseg = makeRandEllipsoid(7, 0.1, 0.1, 0.05, 0.01);
% testseg = makeEllipsoid(9, 0.05, 0.05, 0.05);
% testseg = makeCylinder(9,0.025,0.1);
testseg = makeCone(9,0.1,0.25);

inputstr(1).segobj = testseg;

%% initialice only segment
cmptarget.segments = testseg;
cmptarget.segments.init(inputstr);

%% init inertia
[vol, mas, moi, barycenter, pointcloud, pointmasses] = cmptarget.segments.initSolidInertia(); % O(n^3)
% [vol, mas, moi, pointcloud] = cmptarget.segments.initHollowInertia(); % O(n^3)

%%
close all;
plotSegmentWithPointCloud(cmptarget, pointcloud);
plotSegment(cmptarget);
% plot2tikz('randellipsoid', 1, 'ticks=none');