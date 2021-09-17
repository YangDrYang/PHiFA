function h = plotSegmentWithNormals(cmptarget, pointcloud)
%UNTITLED10 Summary of this function goes here
%   Detailed explanation goes here
h = figure('Name','Target Segment Object Display');
% pcshow(pointcloud,[0 0 0],'MarkerSize',50);
pcshow(pointcloud.');
%             axis tight;
hold on;

% bb = cmptarget.segments.segmentReferenceCube;

% axis( [ bb(2,1) bb(1,1) bb(2,2) bb(1,2) bb(2,3) bb(1,3) ] );

% vbb1 = [bb(1,1) bb(1,2) bb(1,3);...
%     bb(1,1) bb(2,2) bb(1,3);...
%     bb(1,1) bb(2,2) bb(2,3);...
%     bb(1,1) bb(1,2) bb(2,3)];
% f = [ 1 2 3 4 ];
% patch('Faces',f,'Vertices',vbb1,'FaceColor','none','EdgeColor','blue','LineWidth',2);
% vbb2 = [bb(1,1) bb(1,2) bb(1,3);...
%     bb(1,1) bb(1,2) bb(2,3);...
%     bb(2,1) bb(1,2) bb(2,3);...
%     bb(2,1) bb(1,2) bb(1,3)];
% patch('Faces',f,'Vertices',vbb2,'FaceColor','none','EdgeColor','blue','LineWidth',2);
% vbb3 = [bb(1,1) bb(1,2) bb(1,3);...
%     bb(1,1) bb(2,2) bb(1,3);...
%     bb(2,1) bb(2,2) bb(1,3);...
%     bb(2,1) bb(1,2) bb(1,3)];
% patch('Faces',f,'Vertices',vbb3,'FaceColor','none','EdgeColor','blue','LineWidth',2);
% vbb1 = [bb(2,1) bb(1,2) bb(1,3);...
%     bb(2,1) bb(2,2) bb(1,3);...
%     bb(2,1) bb(2,2) bb(2,3);...
%     bb(2,1) bb(1,2) bb(2,3)];
% patch('Faces',f,'Vertices',vbb1,'FaceColor','none','EdgeColor','blue','LineWidth',2);
% vbb2 = [bb(1,1) bb(2,2) bb(1,3);...
%     bb(1,1) bb(2,2) bb(2,3);...
%     bb(2,1) bb(2,2) bb(2,3);...
%     bb(2,1) bb(2,2) bb(1,3)];
% patch('Faces',f,'Vertices',vbb2,'FaceColor','none','EdgeColor','blue','LineWidth',2);
% vbb3 = [bb(1,1) bb(1,2) bb(2,3);...
%     bb(1,1) bb(2,2) bb(2,3);...
%     bb(2,1) bb(2,2) bb(2,3);...
%     bb(2,1) bb(1,2) bb(2,3)];
% patch('Faces',f,'Vertices',vbb3,'FaceColor','none','EdgeColor','blue','LineWidth',2);

f = [1 2 4];
for i = 1:length(cmptarget.segments.facets)
    if size(cmptarget.segments.facets, 2)>1000 && mod(i,100)~=0
            continue;
    end
    vtri = [cmptarget.segments.facets(i).base,...
        cmptarget.segments.facets(i).base+cmptarget.segments.facets(i).edge1,...
        cmptarget.segments.facets(i).base+cmptarget.segments.facets(i).edge1+...
        cmptarget.segments.facets(i).edge2,...
        cmptarget.segments.facets(i).base+cmptarget.segments.facets(i).edge2].';
    patch('Faces',f,'Vertices',vtri,'FaceColor',       [0.8 0.8 1.0], ...
         'EdgeColor',       'black',        ...
         'FaceLighting',    'gouraud',     ...
         'FaceAlpha', 0.1);
end


x(1:length(cmptarget.segments.facets),1) = 0;
y(1:length(cmptarget.segments.facets),1) = 0;
z(1:length(cmptarget.segments.facets),1) = 0;
u(1:length(cmptarget.segments.facets),1) = 0;
v(1:length(cmptarget.segments.facets),1) = 0;
w(1:length(cmptarget.segments.facets),1) = 0;
for i = 1:length(cmptarget.segments.facets)
    x(i) = cmptarget.segments.facets(i).areabarycenter(1);
    y(i) = cmptarget.segments.facets(i).areabarycenter(2);
    z(i) = cmptarget.segments.facets(i).areabarycenter(3);
    u(i) = cmptarget.segments.facets(i).normal(1);
    v(i) = cmptarget.segments.facets(i).normal(2);
    w(i) = cmptarget.segments.facets(i).normal(3);
end

quiver3(x,y,z,u,v,w);

view([120 20]);
hold off;
title('Testobject');
xlabel('x [m]');
ylabel('y [m]');
zlabel('z [m]');
axis equal;
grid;

filename = sprintf('figures\\plot_segment-%s.fig',datestr(now,'yy-mm-dd_HH-MM-SS'));
savefig(filename);
end

