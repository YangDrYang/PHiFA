function out = tophat(x,rRef)
%tophat Summary of this function goes here
%   Detailed explanation goes here
out = zeros(size(x));
for i = 1:numel(x)
    if x(i) < rRef
        out(i) = 1;
    elseif x(i) > rRef
        out(i) = 0;
    else
        out(i) = 0.5;    
    end
end
end

