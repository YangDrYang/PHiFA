function plate = makeSolarPanel(longSide, shortSide, thickness)
%makeCube Creates cube segment with given side length
plate = clCompoundSegment();
for i = 1:6; plate.facets(i) = clRectangle(); end
for i = 1:3; plate.facets(i).base = [0 0 0].'; end
for i = 4:6; plate.facets(i).base = [shortSide longSide thickness].'; end
ex = [1 0 0].'; ey = [0 1 0].'; ez = [0 0 1].';
plate.facets(1).edge1=ez*thickness;
plate.facets(2).edge1=ez*thickness;
plate.facets(3).edge1=ex*shortSide;
plate.facets(4).edge1=-ez*thickness;
plate.facets(5).edge1=-ez*thickness;
plate.facets(6).edge1=-ex*shortSide;

plate.facets(1).edge2=ey*longSide;
plate.facets(2).edge2=ex*shortSide;
plate.facets(3).edge2=ey*longSide;
plate.facets(4).edge2=-ey*longSide;
plate.facets(5).edge2=-ex*shortSide;
plate.facets(6).edge2=-ey*longSide;

plate.facets(1).normal=-ex;
plate.facets(2).normal=-ey;
plate.facets(3).normal=-ez;
plate.facets(4).normal=ex;
plate.facets(5).normal=ey;
plate.facets(6).normal=ez;

plate.facets(1).spflag=true;
plate.facets(2).spflag=false;
plate.facets(3).spflag=true;
plate.facets(4).spflag=true;
plate.facets(5).spflag=false;
plate.facets(6).spflag=true;


for i = 1:6
    plate.facets(i).init();
end

plate.segmentReferenceCube = [longSide shortSide thickness; 0 0 0].';

end