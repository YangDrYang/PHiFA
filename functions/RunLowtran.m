
alt_km = 0;
zenithangle = [0, 30, 85];

p.model=5;
p.h1=alt_km;
p.angle=zenithangle;
p.wlshort= 200;
p.wllong=30000;

T = lowtran_transmission(p);

plot_transmission(T)

hleg = legend([num2str(zenithangle(1)) ' deg.'], [num2str(zenithangle(2)) ' deg.'], ...
    [num2str(zenithangle(3)) ' deg.']);
title(hleg,'Zenith angle');
% hleg.Title.NodeChildren.Position = [0.5 1.5 0];

% plot2tikz('lowtran', 0.7);