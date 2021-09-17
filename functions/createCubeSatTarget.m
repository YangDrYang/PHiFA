function defCube = createCubeSatTarget(units)
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
inputstr.offset = [0; 0; 0];
inputstr.name = 'cube';
inputstr.density = 1000;
inputstr.scale = 1/1000; %units in stl. 1/100 would be centimeter and so on
inputstr.resolution = 250;
inputstr.solid = true;
inputstr.twosided = false;
inputstr.thickness = 0.003;
inputstr.segobj = makePlate(0.1, units*0.1);

defCube = clCompoundTarget();
defCube.name = 'default_CubeSat';
defCube.init(inputstr, tar_input);

end