function seg = makeCone(resolution,radius,height)
%makeCone

seg = clCompoundSegment();
if mod(resolution,2)==0
    resolution=resolution+1;
end

th = linspace(0,2*pi,resolution);
seg.facets(1:2*(resolution-1)) = clTriangle();
fc=1;
for i = 1:resolution-1
    % in sloop
    seg.facets(fc) = clTriangle();
    seg.facets(fc).base = [0 0 height].';
    seg.facets(fc).edge1 = [radius*sin(th(i)) radius*cos(th(i)) 0].' - seg.facets(fc).base;
    seg.facets(fc).edge2 = [radius*sin(th(i+1)) radius*cos(th(i+1)) 0].' - seg.facets(fc).base;
    seg.facets(fc).normal = cross(seg.facets(fc).edge2,...
        seg.facets(fc).edge1);
    seg.facets(fc).init();
    fc=fc+1;
    % bottom
    seg.facets(fc) = clTriangle();
    seg.facets(fc).base = [0 0 0].';
    seg.facets(fc).edge1 = [radius*sin(th(i)) radius*cos(th(i)) 0].' - seg.facets(fc).base;
    seg.facets(fc).edge2 = [radius*sin(th(i+1)) radius*cos(th(i+1)) 0].' - seg.facets(fc).base;
    seg.facets(fc).normal = [0 0 -1].';
    seg.facets(fc).init();
    fc=fc+1;
end

seg.segmentReferenceCube = [radius radius height; -radius -radius 0].';

end

