function seg = makeSolarPanelRight(longSide, shortSide, thickness, resolution)
%makeCube Creates cube segment with given side length
seg = clCompoundSegment();
if mod(resolution,2)==0
    resolution=resolution+1;
end
diameter = 2.8;
seg.segmentReferenceCube = [0 -0.5*diameter 0.5*thickness; -shortSide -0.5*diameter-longSide 0.5*thickness].';
x = zeros(resolution);
y = zeros(resolution);
ori = [-shortSide,-0.5*diameter];
x(end) = 0;
y(end) = -(0.5*diameter+longSide);
for i = 1:resolution-1
    x(i) = ori(1) + shortSide/resolution;
    y(i) = ori(2) + -longSide/resolution;
end

seg.facets(1:2*(resolution-1)) = clRectangle();
fc = 1;

for i = 1:(resolution-1)
    % upper side
    seg.facets(fc) = clTriangle();
    seg.facets(fc).base = [x(i) y(i) 0.5*thickness].';
    seg.facets(fc).edge1 = [x(i+1) y(i) 0.5*thickness].' - seg.facets(fc).base;
    seg.facets(fc).edge2 = [x(i) y(i+1) 0.5*thickness].' - seg.facets(fc).base;
    seg.facets(fc).normal = cross(seg.facets(fc).edge2,...
        seg.facets(fc).edge1);
    fc=fc+1;
    % lower side
    seg.facets(fc) = clTriangle();
    seg.facets(fc).base = [x(i) y(i) -0.5*thickness].';
    seg.facets(fc).edge1 = [x(i+1) y(i) -0.5*thickness].' - seg.facets(fc).base;
    seg.facets(fc).edge2 = [x(i) y(i+1) -0.5*thickness].' - seg.facets(fc).base;
    seg.facets(fc).normal = cross(seg.facets(fc).edge2,...
        seg.facets(fc).edge1);
    fc=fc+1;
end

for i = 1:length(seg.facets)
    seg.facets(i).init();
end
end

