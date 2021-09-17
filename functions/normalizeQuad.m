function norm_quad = normalizeQuad(quad)

norm_quad = zeros(4, size(quad,2));
for i = 1:size(quad,2)
    norm_quad(:,i) = quad(:,i)./norm(quad(:,i));
end
    