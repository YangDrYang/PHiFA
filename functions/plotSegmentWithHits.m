function h = plotSegmentWithHits(cmptarget, hits)
h = plotTarget(cmptarget);

x(1:length(hits),1) = 0;
y(1:length(hits),1) = 0;
z(1:length(hits),1) = 0;
u(1:length(hits),1) = 0;
v(1:length(hits),1) = 0;
w(1:length(hits),1) = 0;
for i = 1:length(hits)
    x(i) = hits(i).hitpos(1);
    y(i) = hits(i).hitpos(2);
    z(i) = hits(i).hitpos(3);
    u(i) = hits(i).impuls(1);
    v(i) = hits(i).impuls(2);
    w(i) = hits(i).impuls(3);
%     u(i) = hits(i).force(1);
%     v(i) = hits(i).force(2);
%     w(i) = hits(i).force(3);
end

quiver3(x,y,z,u,v,w,'LineWidth',2);

view([120 20]);
hold off;
xlabel('x [m]');
ylabel('y [m]');
zlabel('z [m]');
axis equal;
grid;
end

