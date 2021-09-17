function Cn2 = tumodSLCDay(z,za,h0)
%tumodSLCDay Summary of this function goes here
%   Detailed explanation goes here
Cn2(size(z))=0;
for i = 1:numel(z)
    h = getTrueHeight(z(i),za,h0);
    if h < 0
        Cn2(i) = -1;
    elseif h < 18.5
        Cn2(i) = 1.7e-14;
    elseif h < 240
        Cn2(i) = 3.13e-13/h^(1.05);
    elseif h < 880
        Cn2(i) = 1.3e-15;
    elseif h < 7200
        Cn2(i) = 8.87e-7/h^3;
    elseif h < 20000
        Cn2(i) = 2e-16/h^0.5;
    else
        Cn2(i) = 0;
    end
end

end

