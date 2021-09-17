function plate = makePlate(longSide, shortSide)
%makeCube Creates cube segment with given side length
plate = clCompoundSegment();
for i = 1:6; plate.facets(i) = clRectangle(); end
for i = 1:3; plate.facets(i).base = [0 0 0].'; end
for i = 4:6; plate.facets(i).base = [longSide longSide shortSide].'; end
ex = [1 0 0].'; ey = [0 1 0].'; ez = [0 0 1].';
plate.facets(1).edge1=ez*shortSide;
plate.facets(2).edge1=ez*shortSide;
plate.facets(3).edge1=ex*longSide;
plate.facets(4).edge1=-ez*shortSide;
plate.facets(5).edge1=-ez*shortSide;
plate.facets(6).edge1=-ex*longSide;

plate.facets(1).edge2=ey*longSide;
plate.facets(2).edge2=ex*longSide;
plate.facets(3).edge2=ey*longSide;
plate.facets(4).edge2=-ey*longSide;
plate.facets(5).edge2=-ex*longSide;
plate.facets(6).edge2=-ey*longSide;

plate.facets(1).normal=-ex;
plate.facets(2).normal=-ey;
plate.facets(3).normal=-ez;
plate.facets(4).normal=ex;
plate.facets(5).normal=ey;
plate.facets(6).normal=ez;

for i = 1:6
    plate.facets(i).init();
end

plate.segmentReferenceCube = [longSide longSide shortSide; 0 0 0].';

end

