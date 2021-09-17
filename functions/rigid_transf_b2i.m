function z = rigid_transf_b2i(X,EulerAng)
% %±æÃÂµ?ø?º???ªª
% real*8 EulerAng(3) 	%??¿??«
% real*8 X(3),Z(3)	%X±æÃÂ?¯±Í, Zø?º??¯±Í.
% real*8 mattmp1(3),mattmp2(3),rot(3,3)

rot = RotationXYZ(-EulerAng(3),3);
mattmp1 = rot*X;

rot = RotationXYZ(-EulerAng(2),1);
mattmp2=rot*mattmp1;

rot = RotationXYZ(-EulerAng(1),3);
z=rot*mattmp2;

end