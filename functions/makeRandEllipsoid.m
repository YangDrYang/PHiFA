function seg = makeRandEllipsoid(resolution,sidex,sidey,sidez,noise)
%makeEllipsoid Summary of this function goes here
%   Detailed explanation goes here
seg = clCompoundSegment();
if mod(resolution,2)==0
    resolution=resolution+1;
end
[x, y, z] = ellipsoid(0,0,0,sidex/2,sidey/2,sidez/2,resolution);

x1 = x + noise*randn(size(x));
x1(1,:) = x(1,:);
x1(end,:) = x(end,:);
x1(:,1) = x1(:,end);
y1 = y + noise*randn(size(y));
y1(1,:) = y(1,:);
y1(end,:) = y(end,:);
y1(:,1) = y1(:,end);
z1 = z + noise*randn(size(z));
z1(1,:) = z(1,:);
z1(end,:) = z(end,:);
z1(:,1) = z1(:,end);

nfacets = (resolution-1)*2*resolution;
fc = 1; % facetcounter
seg.facets(1:nfacets) = clTriangle();
for i = 1:nfacets; seg.facets(i) = clTriangle(); end

seg.segmentReferenceCube = [max(x1(:)) max(y1(:)) max(z1(:)); min(x1(:)) min(y1(:)) min(z1(:))].';

%first row
for i = 1:resolution
    seg.facets(fc).base = [x1(1,i) y1(1,i) z1(1,i)].';
    seg.facets(fc).edge1 = [x1(2,i) y1(2,i) z1(2,i)].'-seg.facets(fc).base;
    seg.facets(fc).edge2 = [x1(2,i+1) y1(2,i+1) z1(2,i+1)].'-seg.facets(fc).base;
    seg.facets(fc).normal = cross(seg.facets(i).edge2,...
        seg.facets(fc).edge1);
    fc = fc+1;
end

% in between
for j = 2:resolution-1
    for i = 1:resolution
        seg.facets(fc).base = [x1(j,i) y1(j,i) z1(j,i)].';
        seg.facets(fc).edge1 = [x1(j+1,i) y1(j+1,i) z1(j+1,i)].'-seg.facets(fc).base;
        seg.facets(fc).edge2 = [x1(j+1,i+1) y1(j+1,i+1) z1(j+1,i+1)].'-seg.facets(fc).base;
        seg.facets(fc).normal = cross(seg.facets(fc).edge2,...
            seg.facets(fc).edge1);
        fc = fc+1;
        seg.facets(fc).base = [x1(j+1,i+1) y1(j+1,i+1) z1(j+1,i+1)].';
        seg.facets(fc).edge1 = [x1(j,i+1) y1(j,i+1) z1(j,i+1)].'-seg.facets(fc).base;
        seg.facets(fc).edge2 = [x1(j,i) y1(j,i) z1(j,i)].'-seg.facets(fc).base;
        seg.facets(fc).normal = cross(seg.facets(fc).edge2,...
            seg.facets(fc).edge1);
        fc = fc+1;
    end
end

%lastrow
for i = 1:resolution
    seg.facets(fc).base = [x1(end,i) y1(end,i) z1(end,i)].';
    seg.facets(fc).edge1 = [x1(end-1,i) y1(end-1,i) z1(end-1,i)].'-seg.facets(fc).base;
    seg.facets(fc).edge2 = [x1(end-1,i+1) y1(end-1,i+1) z1(end-1,i+1)].'-seg.facets(fc).base;
    seg.facets(fc).normal = cross(seg.facets(fc).edge1,...
        seg.facets(fc).edge2);
    fc = fc+1;
end

for i = 1:nfacets; seg.facets(i).init(); end

end