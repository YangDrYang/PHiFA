wl = 1000E-9;
Ap = 2;
eisert = @(sigma,dist) 1 - exp(-(wl.*dist./Ap+0.5).^2./(8*dist.^2*sigma.^2));
dist = @(pmin, sigma) sqrt(-(2.^2)./(8.*log(1-pmin).*sigma.^2));

dist = [100:100:500 750 1000 1500 2000];
sigma = logspace(-3,3,100)/1000000;

hitprop(1:length(dist),1:length(sigma)) = 0;
dist_plot(1:length(dist),1:length(sigma)) = 0;
sigma_plot(1:length(dist),1:length(sigma)) = 0;

for i = 1:length(dist)
 hitprop(i,:) = eisert(sigma, dist(i)*1000)*100;
 sigma_plot(i,:) = sigma*1000000;
end

for i = 1:length(sigma)
 dist_plot(:, i) = dist;
end

figure;
for i = 1:length(dist)
    semilogx(sigma_plot(i,:), hitprop(i,:),'DisplayName',['dist. ' num2str(dist(i)) 'km']);
    hold on;
end

xlim([5E-2 5E1])
ylim([0 100])
% shading interp
xlabel('Pointing accuracy [murad]');
ylabel('Hit probability [%]');
legend('Location', 'northeastoutside');
grid;
pbaspect([2 1 1])

plot2tikz('hit_prop', 0.66);

%%
% dist = @(pmin, sigma) sqrt(-(2.^2)./(8.*log(1-pmin).*sigma.^2));
% 
% pmin = 0:0.1:1;
% sigma = logspace(-3,3,100)/1000000;
% 
% for i =1:length(pmin)
%     mindist(i,:) = dist(pmin(i),sigma)/1000;
% end
% 
% figure;
% for i = 1:length(pmin)
%     semilogx(sigma.*1E6, mindist(i,:),'DisplayName',['P_min ' num2str(pmin(i)*100) '%']);
%     hold on;
% end
% 
% % xlim([5E-2 5E1])
% % ylim([0 100])
% % shading interp
% xlabel('Pointing accuracy [murad]');
% ylabel('Hit probability [%]');
% legend('Location', 'northeastoutside');
% grid;
% pbaspect([2 1 1])
% 
% plot2tikz('hit_prop', 0.66);
