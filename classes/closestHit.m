function hit = closestHit(hits)
%closestHist function return hit which is closest
hit = hits(1);
    for i = 2:length(hits)
        if hits(i).distFromLaser < hit.distFromLaser
            hit = hits(i);
        end
    end
end

