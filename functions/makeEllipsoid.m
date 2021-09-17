function seg = makeEllipsoid(resolution,sidex,sidey,sidez)
%makeEllipsoid Summary of this function goes here
%   Detailed explanation goes here
seg = clCompoundSegment();
if mod(resolution,2)==0
    resolution=resolution+1;
end
[x, y, z] = ellipsoid(0,0,0,sidex,sidey,sidez,resolution);

nfacets = (resolution-1)*2*resolution;
fc = 1; % facetcounter
seg.facets(1:nfacets) = clTriangle();
for i = 1:nfacets; seg.facets(i) = clTriangle(); end

seg.segmentReferenceCube = [max(x(:)) max(y(:)) max(z(:)); min(x(:)) min(y(:)) min(z(:))].';

%first row
for i = 1:resolution
    seg.facets(fc).base = [x(1,i) y(1,i) z(1,i)].';
    seg.facets(fc).edge1 = [x(2,i) y(2,i) z(2,i)].'-seg.facets(fc).base;
    seg.facets(fc).edge2 = [x(2,i+1) y(2,i+1) z(2,i+1)].'-seg.facets(fc).base;
    seg.facets(fc).normal = cross(seg.facets(i).edge2,...
        seg.facets(fc).edge1);
    fc = fc+1;
end

% in between
for j = 2:resolution-1
    for i = 1:resolution
        seg.facets(fc).base = [x(j,i) y(j,i) z(j,i)].';
        seg.facets(fc).edge1 = [x(j+1,i) y(j+1,i) z(j+1,i)].'-seg.facets(fc).base;
        seg.facets(fc).edge2 = [x(j+1,i+1) y(j+1,i+1) z(j+1,i+1)].'-seg.facets(fc).base;
        seg.facets(fc).normal = cross(seg.facets(fc).edge2,...
            seg.facets(fc).edge1);
        fc = fc+1;
        seg.facets(fc).base = [x(j+1,i+1) y(j+1,i+1) z(j+1,i+1).'];
        seg.facets(fc).edge1 = [x(j,i+1) y(j,i+1) z(j,i+1)].'-seg.facets(fc).base;
        seg.facets(fc).edge2 = [x(j,i) y(j,i) z(j,i)].'-seg.facets(fc).base;
        seg.facets(fc).normal = cross(seg.facets(fc).edge2,...
            seg.facets(fc).edge1);
        fc = fc+1;
    end
end

%lastrow
for i = 1:resolution
    seg.facets(fc).base = [x(end,i) y(end,i) z(end,i)].';
    seg.facets(fc).edge1 = [x(end-1,i) y(end-1,i) z(end-1,i)].'-seg.facets(fc).base;
    seg.facets(fc).edge2 = [x(end-1,i+1) y(end-1,i+1) z(end-1,i+1)].'-seg.facets(fc).base;
    seg.facets(fc).normal = cross(seg.facets(fc).edge1,...
        seg.facets(fc).edge2);
    fc = fc+1;
end

for i = 1:nfacets; seg.facets(i).init(); end
end

