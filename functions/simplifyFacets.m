function outfac = simplifyFacets(infac, maxang, maxdist)
nOutFacCount = 0;
outfac(1:length(infac)) = clRectangle();
alreadyAllocated(1:length(infac)) = 0;
for i = 1:length(infac)
    for j = 1:length(infac)
        if i == j || ~isempty(find(alreadyAllocated==i)) || ~isempty(find(alreadyAllocated==j))
            continue;
        end
        edge(1:3,1:3) = 0;
        base = [0;0;0];
        if acos(dot(infac(i).normal, infac(j).normal)) < maxang
%             fprintf('base-base %f\n', norm(infac(i).base-infac(j).base));
%             fprintf('base-edge1 %f\n', norm(infac(i).base-infac(j).base-infac(j).edge1));
%             fprintf('base-edge2 %f\n', norm(infac(i).base-infac(j).base-infac(j).edge2));
%             fprintf('edge1-base %f\n', norm(infac(i).base+infac(i).edge1-infac(j).base));
%             fprintf('edge1-edge1 %f\n', norm(infac(i).base+infac(i).edge1-infac(j).base-infac(j).edge1));
%             fprintf('edge1-edge2 %f\n', norm(infac(i).base+infac(i).edge1-infac(j).base-infac(j).edge2));
%             fprintf('edge2-base %f\n', norm(infac(i).base+infac(i).edge2-infac(j).base));
%             fprintf('edge2-edge1 %f\n', norm(infac(i).base+infac(i).edge2-infac(j).base-infac(j).edge1));
%             fprintf('edge2-edge2 %f\n', norm(infac(i).base+infac(i).edge2-infac(j).base-infac(j).edge2));
%             fprintf('base1 %s\n', sprintf('%f ', infac(i).base));
%             fprintf('edge11 %s\n',sprintf('%f ', infac(i).base+infac(i).edge1));
%             fprintf('edge21 %s\n', sprintf('%f ', infac(i).base+infac(i).edge2));
%             fprintf('base2 %s\n', sprintf('%f ', infac(j).base));
%             fprintf('edge12 %s\n', sprintf('%f ', infac(j).base+infac(j).edge1));
%             fprintf('edge22 %s\n', sprintf('%f ', infac(j).base+infac(j).edge2));
            
            sameEdgesCount = 0;
            if norm(infac(i).base-infac(j).base)<maxdist || ...
                    norm(infac(i).base-infac(j).base-infac(j).edge1)<maxdist || ...
                    norm(infac(i).base-infac(j).base-infac(j).edge2)<maxdist
                sameEdgesCount = sameEdgesCount+1;
                edge(:,sameEdgesCount) = infac(i).base;
            else
                base = infac(i).base;
            end
            
            if norm(infac(i).base+infac(i).edge1-infac(j).base)<maxdist || ...
                    norm(infac(i).base+infac(i).edge1-infac(j).base-infac(j).edge1)<maxdist || ...
                    norm(infac(i).base+infac(i).edge1-infac(j).base-infac(j).edge2)<maxdist
                sameEdgesCount = sameEdgesCount+1;
                edge(:,sameEdgesCount) = infac(i).base+infac(i).edge1;
            else
                base = infac(i).base+infac(i).edge1;
            end
            
            if norm(infac(i).base+infac(i).edge2-infac(j).base)<maxdist || ...
                    norm(infac(i).base+infac(i).edge2-infac(j).base-infac(j).edge1)<maxdist || ...
                    norm(infac(i).base+infac(i).edge2-infac(j).base-infac(j).edge2)<maxdist
                sameEdgesCount = sameEdgesCount+1;
                edge(:,sameEdgesCount) = infac(i).base+infac(i).edge2;
            else
                base = infac(i).base+infac(i).edge2;
            end
            
            edge(:,1) = edge(:,1) - base;
            edge(:,2) = edge(:,2) - base;
            antibase = base + edge(:,1) + edge(:,2);
%             fprintf('base-base %f\n', norm(antibase-infac(j).base));
%             fprintf('base-edge1 %f\n', norm(antibase-infac(j).base-infac(j).edge1));
%             fprintf('base-edge2 %f\n', norm(antibase-infac(j).base-infac(j).edge2));
            if sameEdgesCount == 2 && ...
                    ( norm(antibase-infac(j).base)<maxdist || ...
                    norm(antibase-infac(j).base-infac(j).edge1)<maxdist || ...
                    norm(antibase-infac(j).base-infac(j).edge2)<maxdist )
                nOutFacCount = nOutFacCount + 1;
                outfac(nOutFacCount) = clRectangle();
                outfac(nOutFacCount).base = base;
                outfac(nOutFacCount).edge1 = edge(:,1);
                outfac(nOutFacCount).edge2 = edge(:,2);
                outfac(nOutFacCount).normal = infac(i).normal;
                outfac(nOutFacCount).init();
                alreadyAllocated(2*nOutFacCount-1) = i;
                alreadyAllocated(2*nOutFacCount) = j;
            end
        end
    end
end

outfac = outfac(1:nOutFacCount);

end