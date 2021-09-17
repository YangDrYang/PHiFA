function bb = getBoundingBoxfromFacets(facets)
bb = [-1000 -1000 -1000; 1000 1000 1000].';

for i=1:length(facets)
    bb = facets(i).compareReferenceCube(bb);
end

end