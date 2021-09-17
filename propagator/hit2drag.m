function drag = hit2drag(hit)

drag = clDrag();

drag.hitpos = hit.hitpos;
drag.las2hitpos = hit.hitpos;
drag.distFromLaser = hit.hitpos;

end