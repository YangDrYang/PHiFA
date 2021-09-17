function seg = makeTruncatedCone(resolution, radiusFront, radiusEnd, len)
%makeTruncatedCone 

seg = clCompoundSegment();
if mod(resolution,2)==0
    resolution=resolution+1;
end

th = linspace(0,2*pi,resolution);
seg.facets(1:(resolution-1)) = clRectangle();
fc=1;
for i = 1:resolution-1
    % in sloop
    seg.facets(fc) = clRectangle();
    seg.facets(fc).base = [0 radiusEnd*sin(th(i)) radiusEnd*cos(th(i))].';
    seg.facets(fc).edge1 = [0 radiusEnd*sin(th(i+1)) radiusEnd*cos(th(i+1))].' - seg.facets(fc).base;
    seg.facets(fc).edge2 = [len, radiusFront*sin(th(i)) radiusFront*cos(th(i+1))].' - seg.facets(fc).base;
    seg.facets(fc).normal = cross(seg.facets(fc).edge2,...
        seg.facets(fc).edge1);
    seg.facets(fc).init();
    fc=fc+1;
end

seg.segmentReferenceCube = [len max(radiusFront,radiusEnd) max(radiusFront,radiusEnd); 0 -max(radiusFront,radiusEnd) -max(radiusFront,radiusEnd)].';

end

