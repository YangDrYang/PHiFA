    function [intmax,poserror] = findMaxIntensity(Uin,x,y)
%findMaxIntensity Finds position of the maximum intensity of given
%distribution
Uinabs = abs(Uin);
[intmax, posmax] = max(Uinabs(:));
[xmax,ymax] = ind2sub(size(Uin),posmax);
poserror = sqrt(x(xmax,ymax)^2 + y(xmax,ymax)^2);
end

