function [facets,nfacets] = STL2Facets(stlfile)
%STL2Facets Summary of this function goes here
%   Detailed explanation goes here
            disp('Reading STL input file');
            [f,v,n] = stlread(stlfile);
%             obj.facets(size(f,1)) = clFacet(true, eFacetShape.plate);
            disp('Transforming to facets');
            nfacets = size(f,1);
            facets(1,nfacets) = clTriangle();
            for i = 1:size(f,1)
                facets(i) = clTriangle();
                facets(i).base = v(f(i,1), :)/1000;
                facets(i).edge1 = (v(f(i,2), :) - v(f(i,1), :))/1000;
                facets(i).edge2 = (v(f(i,3), :) - v(f(i,1), :))/1000;
                facets(i).normal = n(i,:)/abs(n(i,:));
            end
end

