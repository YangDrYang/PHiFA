function qout = norm_quat(qin)
    norm = sqrt(qin(1).^2+qin(2).^2+qin(3).^2+qin(4).^2);
    qout = qin./norm;
end