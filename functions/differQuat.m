function dq = differQuat(q, w)

Omega = [0, -w(1), -w(2), -w(3); ...
    w(1), 0, w(3), -w(2); ...
    w(2), -w(3), 0, w(1); ...
    w(3), w(2), -w(1), 0];

dq = 0.5*( Omega*q );

end