function cube = makeCube(sideLength)
%makeCube Creates cube segment with given side length
cube = clCompoundSegment();
for i = 1:6; cube.facets(i) = clRectangle(); end
for i = 1:3; cube.facets(i).base = [-sideLength/2 -sideLength/2 0].'; end
for i = 4:6; cube.facets(i).base = [sideLength/2 sideLength/2 sideLength].'; end
ex = [1 0 0].'; ey = [0 1 0].'; ez = [0 0 1].';
cube.facets(1).edge1=ez*sideLength;
cube.facets(2).edge1=ez*sideLength;
cube.facets(3).edge1=ex*sideLength;
cube.facets(4).edge1=-ez*sideLength;
cube.facets(5).edge1=-ez*sideLength;
cube.facets(6).edge1=-ex*sideLength;

cube.facets(1).edge2=ey*sideLength;
cube.facets(2).edge2=ex*sideLength;
cube.facets(3).edge2=ey*sideLength;
cube.facets(4).edge2=-ey*sideLength;
cube.facets(5).edge2=-ex*sideLength;
cube.facets(6).edge2=-ey*sideLength;

cube.facets(1).normal=-ex;
cube.facets(2).normal=-ey;
cube.facets(3).normal=-ez;
cube.facets(4).normal=ex;
cube.facets(5).normal=ey;
cube.facets(6).normal=ez;

for i = 1:6
    cube.facets(i).init();
end

cube.segmentReferenceCube = [sideLength/2 sideLength/2 sideLength; -sideLength/2 -sideLength/2 0].';

end

