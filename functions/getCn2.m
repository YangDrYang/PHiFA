function Cn2 = getCn2(tumod,za,h0,dist)
%intCn2dz Function integrates Cn2 times delta z over given distance

Cn2(size(dist))=0;
if dist(1)~=0
    dist = [0 dist];
else
    Cn2 = Cn2(1:end-1);
end

n=length(dist);
for i = 1:n-1
    Cn2(i) = integral(@(x)tumod(x,za,h0),dist(i),dist(i+1),'ArrayValued',true)/(dist(i+1)-dist(i));
end

end

