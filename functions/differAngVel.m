function dw = differAngVel(moi, moi_inv, w, g)

Iw = moi*w;

cwIwg = [ -w(2)*Iw(3)+w(3)*Iw(2)+g(1); ...
        w(1)*Iw(3)-w(3)*Iw(1)+g(2); ...
        -w(1)*Iw(2)+w(2)*Iw(1)+g(3) ];
    
dw = moi_inv*cwIwg;

end