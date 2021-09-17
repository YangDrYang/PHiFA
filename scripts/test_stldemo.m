%% 3D Model Demo
% This is short demo that loads and renders a 3D model of a human femur. It
% showcases some of MATLAB's advanced graphics features, including lighting and
% specular reflectance.

% Copyright 2011 The MathWorks, Inc.


%% Load STL mesh
% Stereolithography (STL) files are a common format for storing mesh data. STL
% meshes are simply a collection of triangular faces. This type of model is very
% suitable for use with MATLAB's PATCH graphics object.

% Import an STL mesh, returning a PATCH-compatible face-vertex structure
% [f,v,n] = stlread('femur.stl');
% [f,v,n] = stlread('67P-quad.stl');
% [f,v,n] = stlread('COMET_67P_C-G.stl');
% [f,v,n] = stlread('stlfiles\femur.stl');
% stlfile = 'stlfiles/CubeSat_middle.STL';
% stlfile = 'stlfiles/3UCubeSat.stl';
stlfile = 'stlfiles/3UCubeSat2.stl';
[f,v,n] = stlread(stlfile);

%% Render
% The model is rendered with a PATCH graphics object. We also add some dynamic
% lighting, and adjust the material properties to change the specular
% highlighting.
fv = struct('faces',f,'vertices',v);
patch(fv,'FaceColor',       [0.8 0.8 1.0], ...
         'EdgeColor',       'none',        ...
         'FaceLighting',    'gouraud',     ...
         'AmbientStrength', 0.15,          ...
         'FaceAlpha', 0.1);

% Add a camera light, and tone down the specular highlighting
camlight('headlight');
material('dull');

% Fix the axes scaling, and set a nice view angle
axis('image');
% view([-135 35]);
view([120 20]);
plot2tikz('3ucubesat_stl',0.6);

[facets,~] = STL2Facets(stlfile);
segment.facets = facets;
plotSegment(segment);