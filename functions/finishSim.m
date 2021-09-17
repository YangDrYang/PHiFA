propagator.finishSim();
save2file = sprintf('%s.mat', propagator.outfilename);
if exist('eph', 'var')
    save(save2file, 'propagator', 'initstate', 'eph', 'Y0', 'timestamp');
else
    save(save2file, 'propagator', 'initstate', 'Y0', 'timestamp');
end