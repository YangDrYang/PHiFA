function seg = makeCylinderFront(resolution, radius, height)
%makeCylinder 
seg = clCompoundSegment();
if mod(resolution,2)==0
    resolution=resolution+1;
end
[x,y,z] = cylinder(radius, resolution-1);
x0 = z;
y0 = -x;
z0 = -y;
x = x0;
y = y0;
z = z0;

x = x.*height;
seg.segmentReferenceCube = [max(x(:)) max(y(:)) max(z(:)); min(x(:)) min(y(:)) min(z(:))].';

seg.facets(1:2*(resolution-1)) = clRectangle();
fc = 1;

for i = 1:(resolution-1)
%     % bottom
%     seg.facets(fc) = clTriangle();
%     seg.facets(fc).base = [x(1,i) 0 0].';
%     seg.facets(fc).edge1 = [x(1,i) y(1,i) z(1,i)].' - seg.facets(fc).base;
%     seg.facets(fc).edge2 = [x(1,i+1) y(1,i+1) z(1,i+1)].' - seg.facets(fc).base;
%     seg.facets(fc).normal = cross(seg.facets(fc).edge2,...
%         seg.facets(fc).edge1);
%     fc=fc+1;
    % in between
    seg.facets(fc) = clRectangle();
    seg.facets(fc).base = [x(1,i) y(1,i) z(1,i)].';
    seg.facets(fc).edge1 = [x(2,i) y(2,i) z(2,i)].' - seg.facets(fc).base;
    seg.facets(fc).edge2 = [x(1,i+1) y(1,i+1) z(1,i+1)].' - seg.facets(fc).base;
    seg.facets(fc).normal = cross(seg.facets(fc).edge2,...
        seg.facets(fc).edge1);
    fc=fc+1;
    % top
    seg.facets(fc) = clTriangle();
    seg.facets(fc).base = [x(2,i) 0 0].';
    seg.facets(fc).edge1 = [x(2,i) y(2,i) z(2,i)].' - seg.facets(fc).base;
    seg.facets(fc).edge2 = [x(2,i+1) y(2,i+1) z(2,i+1)].' - seg.facets(fc).base;
    seg.facets(fc).normal = cross(seg.facets(fc).edge1,...
        seg.facets(fc).edge2);
    fc=fc+1;
end

for i = 1:length(seg.facets)
    seg.facets(i).init();
end

end

