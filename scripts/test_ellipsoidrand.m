[x, y, z] = ellipsoid(0,0,0,8,15,5,5);

fvc = surf2patch(x,y,z,z);

%% plot
figure;
x1 = x + rand(size(x));
x1(1,:) = x(1,:);
x1(end,:) = x(end,:);
x1(:,1) = x1(:,end);
% x1(:,1) = x(:,1);
% x1(:,end) = x(:,end);
y1 = y + rand(size(y));
y1(1,:) = y(1,:);
y1(end,:) = y(end,:);
y1(:,1) = y1(:,end);
% y1(:,1) = y(:,1);
% y1(:,end) = y(:,end);
z1 = z + rand(size(z));
z1(1,:) = z(1,:);
z1(end,:) = z(end,:);
z1(:,1) = z1(:,end);

surf(x1, y1, z1)
axis equal

figure;
surf(x, y, z)
axis equal