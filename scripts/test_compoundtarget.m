%test script

clear;
%% initialice shape

inputstr = clCompoundSegmentsInput;

% inputstr(1).bSTL = true;
% inputstr(1).stlFile = string('stlfiles\tetrahedron.stl');
% inputstr.stlFile = string('stlfiles\femur.stl');

inputstr(1).bSTL = false;
inputstr(1).stlFile = "";

inputstr(1).surfAtr = clSA_wilken2015();
inputstr(1).surfAtr.init([1 6 3 4 5 6], [1 2 3 4]);
inputstr(1).surfAtr.name = 'testsurf1';
inputstr(1).offset = [0; 0; 0];
inputstr(1).name = 'test1';
inputstr(1).density = 2500;
inputstr(1).scale = 1/1000; %units in stl. 1/100 would be centimeter and so on
inputstr(1).resolution = 300;
inputstr(1).solid = true;
inputstr(1).twosided = false;
inputstr(1).thickness = 0.001;

inputstr(1).segobj = makeCylinder(9,0.025,0.1);

% inputstr(2).bSTL = true;
% inputstr(2).stlFile = string('stlfiles\tetrahedron.stl');
% inputstr(2).stlFile = string('stlfiles\femur.stl');

inputstr(2).bSTL = false;
inputstr(2).stlFile = "";

inputstr(2).surfAtr = clSA_wilken2015();
inputstr(2).surfAtr.init([1 2 3 4 5 6], [1 2 3 4]);
inputstr(2).surfAtr.name = 'testsurf2';
inputstr(2).offset = [0; 0; 0.1];
inputstr(2).name = 'test2';
inputstr(2).density = 2700;
inputstr(2).scale = 1/1000; %units in stl. 1/100 would be centimeter and so on
inputstr(2).resolution = 300;
inputstr(2).solid = true;
inputstr(2).twosided = false;
inputstr(2).thickness = 0.001;

inputstr(2).segobj = makeCone(9,0.025,0.07);

% inputstr(2).bSTL = true;
% inputstr(2).stlFile = string('stlfiles\tetrahedron.stl');
% inputstr(2).stlFile = string('stlfiles\femur.stl');

inputstr(3).bSTL = false;
inputstr(3).stlFile = "";

inputstr(3).surfAtr = clSA_wilken2015();
inputstr(3).surfAtr.init([1 2 3 4 5 6], [1 2 3 4]);
inputstr(3).surfAtr.name = 'testsurf2';
inputstr(3).offset = [0; 0; -0.05];
inputstr(3).name = 'test2';
inputstr(3).density = 2700;
inputstr(3).scale = 1/1000; %units in stl. 1/100 would be centimeter and so on
inputstr(3).resolution = 300;
inputstr(3).solid = true;
inputstr(3).twosided = false;
inputstr(3).thickness = 0.001;

inputstr(3).segobj = makeCube(0.05);

tmpinput = clCompoundTargetInput();

%% initialice target
cmptarget = clCompoundTarget();
cmptarget.init(inputstr, tmpinput);

%% initialice only segment
% cmptarget.segments.init(inputstr);
% [vol, mas, moi, pointcloud] = cmptarget.segments.initInertia(); % O(n^3)

%% Output
plotTarget(cmptarget);

%%
% plot2tikz('cylinderconecube', 1, 'ticks=none');
