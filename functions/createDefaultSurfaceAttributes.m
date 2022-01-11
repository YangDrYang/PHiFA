function sa = createDefaultSurfaceAttributes(tau,wl,A)
    sa = clSA_lorbeer2018experimental();
    sa.init(tau,wl,A);
    sa.name = 'aluminium';
end