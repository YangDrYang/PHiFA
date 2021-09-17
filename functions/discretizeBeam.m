function [rays, spacing] = discretizeBeam(ray, radius, resolution)
%discretizeBeam Function takes ray, area and resolution and returns number
%of rays parallel to each other covering area with given resolution.
%Might be improved by projecting bounding box on plane orthogonal to ray
%direction

% ray defined in rigid body coordinate system

l_origin = norm(ray.origin);
nelements = ceil(2*radius*resolution);

if mod(nelements,2)==0
    nelements = nelements + 1;
end
increment = 1/resolution;
maxradius=radius+1/resolution;

rays(1:nelements^2) = clRay();
nrays = 0;
stpt = [-nelements/2*increment -nelements/2*increment];
for i = 1:nelements
    for j = 1:nelements
        ipt = stpt + [i*increment j*increment];
        if norm(ipt)<maxradius
            nrays=nrays+1;
            rays(nrays) = clRay();
            rays(nrays).direction = ray.direction;
            rays(nrays).origin = ray.xdir*ipt(1) + ray.ydir*ipt(2) - l_origin*ray.direction;
            rays(nrays).xdir = ray.xdir; rays(nrays).ydir = ray.ydir;
        end
    end
end
rays = rays(1:nrays);
spacing = increment;
end

