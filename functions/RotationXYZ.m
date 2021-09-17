%Rotation Matrix ????æÿ?Û£®XYZ£©

function Rot = RotationXYZ(theta,Nxyz)
% real*8 theta,Rot(3,3),cost,sint
% integer Nxyz	!1£¨2£¨3?÷±????¶»?X£¨Y£¨Z÷·????

cost=cos(theta);
sint=sin(theta);

Rot=zeros(3,3);

switch Nxyz
    case 1
        Rot(1,1)=1;
        Rot(2,2)=cost;
        Rot(2,3)=sint;
        Rot(3,2)=-sint;
        Rot(3,3)=cost;
    case 2
        Rot(1,1)=cost;
        Rot(1,3)=-sint;
        Rot(2,2)=1;
        Rot(3,1)=sint;
        Rot(3,3)=cost;
    case 3
        Rot(1,1)=cost;
        Rot(1,2)=sint;
        Rot(2,1)=-sint;
        Rot(2,2)=cost;
        Rot(3,3)=1;
end
end