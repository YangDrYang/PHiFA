function spinVal = spinPara(angVel)

%%from angular velocity to spin parameters: magnitude of angular velocity,
%%azimuth and elevation angles

omega = zeros(1,size(angVel,2));
Az = zeros(1,size(angVel,2));
El = zeros(1,size(angVel,2));
for i = 1:size(angVel,2)
    omega(i) = norm(angVel(1:3,i));
    Az(i) = asin(angVel(3,i)/omega(i));
    El(i) = atan2(angVel(2,i),angVel(1,i));
end
spinVal = [omega;Az;El];