function M = TRIC(r,v)
%r 3*1 vector
%v 3*1 vector
u3 = cross(r,v)./norm(cross(r,v));%3*1 vector
u1 = r./norm(r);%3*1 vector
u2 = cross(u3,u1);%3*1 vector

M = [u1,u2,u3];