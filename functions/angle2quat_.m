function q = angle2quat(ang)

cang(1:3)=0;
for i = 1:3
    cang(i) = cos(ang(i)/2.0);
end
sang(1:3)=0;
for i = 1:3
    sang(i) = sin(ang(i)/2.0); 
end

q(4,1) = 0;
q(1) =     cang(1)*cang(2)*cang(3)+sang(1)*sang(2)*sang(3);
q(2) =     cang(1)*cang(2)*sang(3)-sang(1)*sang(2)*cang(3);
q(3) =     cang(1)*sang(2)*cang(3)+sang(1)*cang(2)*sang(3);
q(4) =     sang(1)*cang(2)*cang(3)-cang(1)*sang(2)*sang(3);

end