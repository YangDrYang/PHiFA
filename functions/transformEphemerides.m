function eph_out = transformEphemerides(trans_fcn, eph, start_mjd)
    eph_out = zeros(size(eph));
    eph_out(1,:) = eph(1,:);
    for i = 1:length(eph)
        eph_out(2:7,i) = trans_fcn(start_mjd+eph_out(1,i)/86400, eph(2:7,i)')';
    end
end