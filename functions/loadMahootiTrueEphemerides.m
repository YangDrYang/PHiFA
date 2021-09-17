function true_eph = loadMahootiTrueEphemerides()
True_EnvisatStates;
True_Eph = True_Eph';
true_eph = zeros(size(True_Eph)+[1 0]);
true_eph(2:end,:) = True_Eph;
t = 0;
for i = 1:length(true_eph)
    true_eph(1,i) = t;
    t=t+60;
end
end