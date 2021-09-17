function propagator = reconstructPerturbations(propagator, ephemerides)

propagator = clPropagator.instance(propagator);

for i = 1:length(ephemerides)
    propagator.i_step=i;
    propagator.Accel(ephemerides(1,i),ephemerides(2:14,i));
end

end