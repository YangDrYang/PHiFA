% sc_plotCoCo

wvls = 400:200:1400;
wvls = wvls.*1E-9;
pl = 1E-9;

pls = [0.1, 0.25, 0.5, 1, 2.5, 10];
pls = pls.*1E-9;
wl = 1064E-9;

A = 26.98;

fluence = logspace(-1, 2);
fluence = fluence.*1E4;

cocowvls = zeros(length(fluence), length(wvls));
cocopls = zeros(length(fluence), length(pls));

for i = 1:6
    sawls(i) = createDefaultSurfaceAttributes(pl, wvls(i), A);
    sawls(i).bUseDeviation = true;
    sapls(i) = createDefaultSurfaceAttributes(pls(i), wl, A);
    sapls(i).bUseDeviation = true;
end

for i = 1:6
    for j = 1:length(fluence)
        cocowvls(j,i) = sawls(i).getCouplingCoeffientAblation(fluence(j));
        cocopls(j,i) = sapls(i).getCouplingCoeffientAblation(fluence(j));
    end
end

fluence = fluence.*1E-4;
cocowvls = cocowvls.*1E6;
cocopls = cocopls.*1E6;
figure;
for i = 1:6
    semilogx(fluence(:), cocowvls(:,i),'DisplayName',[num2str(wvls(i)*1E9) ' nm']);
    hold on;
end
xlabel('Fluence [J/cm^2]','FontSize',14);
ylabel('Coupling Coefficient [\muN/W]','FontSize',14);
xlim([1E-1 1E2]);
ylim([0 35]);
legend('Location', 'northwest','FontSize',14);
grid;

% plot2tikz('cocowvls',0.8);

figure;
for i = 1:6
    semilogx(fluence(:), cocopls(:,i), 'DisplayName',[num2str(pls(i)*1E9) ' ns']);
    hold on;
end
xlabel('Fluence [J/cm^2]','FontSize',14);
ylabel('Coupling Coefficient [\muN/W]','FontSize',14);
ylim([0, 35]);
legend('Location', 'northeast','FontSize',14);
grid;

% plot2tikz('cocopls',0.8);