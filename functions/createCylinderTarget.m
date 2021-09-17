function defCube = createCylinderTarget(height, bSolid)
%createDefaultCube Function creates target comprised of one hollow cube
%with edge length of 10 cm similar to cubesat
tar_input = clCompoundTargetInput();

inputstr.bSTL = false;
inputstr.stlFile = "";

inputstr.surfAtr = clSA_wilken2015();
a = [-0.01314 0.01753 0.01904 449.81 632.06 913.51];
t = [15.04 1.19 2.41 23.38];
inputstr.surfAtr.init(a, t);
inputstr.surfAtr.name = 'wilken2015_default';
inputstr.offset = [0; 0; -5.0];
inputstr.name = 'default_Cube';
inputstr.density = 2700;
inputstr.scale = 1/1000; %units in stl. 1/100 would be centimeter and so on
inputstr.resolution = 25/height;
inputstr.solid = bSolid;
inputstr.twosided = false;
inputstr.thickness = 0.003;
inputstr.segobj = makeCylinder(20,height/4,height);

% nSegments = length(inputstr);
defCube = clCompoundTarget();
defCube.name = 'Test Cylinder';
% defCube.init(nSegments, inputstr);
defCube.init(inputstr,tar_input);

end